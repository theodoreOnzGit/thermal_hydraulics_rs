use std::{f64::consts::PI, thread::JoinHandle};
use std::thread;

use uom::ConstZero;
use uom::si::pressure::atmosphere;
use uom::si::f64::*;
use ndarray::*;
use super::InsulatedPipe;
use crate::boussinesq_solver::heat_transfer_correlations::thermal_resistance::try_get_thermal_conductance_annular_cylinder;
use crate::boussinesq_solver::{heat_transfer_correlations::nusselt_number_correlations::{enums::NusseltCorrelation, input_structs::GnielinskiData}, pre_built_components::heat_transfer_entities::preprocessing::try_get_thermal_conductance_based_on_interaction};
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boundary_conditions::BCType;
use crate::boussinesq_solver::pre_built_components::heat_transfer_entities::HeatTransferEntity;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidComponentTrait;

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

impl InsulatedPipe {


    /// used to connect the arrays laterally 
    /// you'll need to set the mass flowrate and heater power
    ///
    /// executes serially, and uses lots of cloning, so it's 
    /// heavier in resource usage,
    ///
    /// unoptimised in this regard
    /// at each timestep, you are allowed to set a heater power, where 
    /// heat is dumped into the heated tube surrounding the pipe
    /// you set it using the heater power input here.
    ///
    /// otherwise you set it to zero for an unpowered pipe
    #[inline]
    pub fn lateral_and_miscellaneous_connections(&mut self,
        mass_flowrate: MassRate,
        heater_power: Power) -> Result<(), ThermalHydraulicsLibError>{


        // first let's get all the conductances 
        
        // |                        |               |               |
        // |                        |               |               |
        // |------ fluid -----------|----shell------|----insulation-| ambient
        // |                        |               |               |
        // |                        |               |               |
        //
        // 1. we'll need the ambient to insulation midpoint (nodal) thermal conductance
        let heat_transfer_to_ambient: HeatTransfer = self.heat_transfer_to_ambient;

        let insulation_to_air_nodal_conductance: ThermalConductance 
        = self.get_ambient_surroundings_to_insulation_thermal_conductance(
            heat_transfer_to_ambient
        )?;

        // 2. we'll need the fluid to shell midpoint thermal conductance

        self.set_mass_flowrate(mass_flowrate);

        let pipe_shell_node_to_fluid_array_conductance: ThermalConductance 
        = self.get_fluid_array_to_pipe_shell_conductance_no_wall_temp_correction()?;

        // 3. we'll need the shell midpoint to insulation midpoint thermal conductance

        let pipe_shell_to_insulation_array_conductance: ThermalConductance 
            = self.get_pipe_shell_to_insulation_conductance()?;

        // next, we need to consider discretisation, ie how much 
        // power fraction
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let q_fraction_per_node: f64 = 1.0/ number_of_temperature_nodes as f64;
        let mut q_frac_arr: Array1<f64> = Array::default(number_of_temperature_nodes);
        q_frac_arr.fill(q_fraction_per_node);

        // then get the ambient temperature 

        let ambient_temp = self.ambient_temperature;

        // lateral connections 
        {
            // first i will need to create temperature vectors 
            // for ambient temperature. This is for use in calculating 
            // heat loss from insulation to ambient air

            let mut ambient_temperature_vector: Vec<ThermodynamicTemperature> 
            = Array1::default(number_of_temperature_nodes)
                .iter().map( |&temp| {
                    temp
                }
                ).collect();

            ambient_temperature_vector.fill(ambient_temp);


            // for this process, I will make a clone of 
            // each HeatTransferEntity, modify the clone, then 
            // replace the HeatTransferEntity within the pipe using 
            // these changed entities

            let mut pipe_shell_clone: SolidColumn = 
                self.pipe_shell.clone().try_into()?;

            let mut fluid_array_clone: FluidArray = 
                self.pipe_fluid_array.clone().try_into()?;

            let mut insulation_array_clone: SolidColumn = 
                self.insulation.clone().try_into()?;



            // note, must set mass flowrate first 
            // otherwise there is by default zero flow through 
            // the array

            fluid_array_clone.set_mass_flowrate(
                mass_flowrate);

            // temperature vectors

            let pipe_temp_vector: Vec<ThermodynamicTemperature> 
            = pipe_shell_clone.get_temperature_vector()?;

            let fluid_temp_vector: Vec<ThermodynamicTemperature> 
            = fluid_array_clone.get_temperature_vector()?;

            let insulation_temp_vector: Vec<ThermodynamicTemperature> 
            = insulation_array_clone.get_temperature_vector()?;

            // second, fill them into the each array 
            
            // insulation to air interaction

            insulation_array_clone.lateral_link_new_temperature_vector_avg_conductance(
                insulation_to_air_nodal_conductance,
                ambient_temperature_vector
            )?;

            // insulation to shell interaction 

            pipe_shell_clone.lateral_link_new_temperature_vector_avg_conductance(
                pipe_shell_to_insulation_array_conductance,
                insulation_temp_vector.clone()
            )?;

            insulation_array_clone.lateral_link_new_temperature_vector_avg_conductance(
                pipe_shell_to_insulation_array_conductance,
                pipe_temp_vector.clone()
            )?;

            // pipe_shell shell to fluid_array interaction

            pipe_shell_clone.lateral_link_new_temperature_vector_avg_conductance(
                pipe_shell_node_to_fluid_array_conductance,
                fluid_temp_vector.clone()
            )?;

            fluid_array_clone.lateral_link_new_temperature_vector_avg_conductance(
                pipe_shell_node_to_fluid_array_conductance,
                pipe_temp_vector
            )?;

            // we also want to add a heat source to pipe_shell shell
            // this is zero by default, but I'm leaving it here as a 
            // scaffold for heated pipes

            pipe_shell_clone.lateral_link_new_power_vector(
                heater_power,
                q_frac_arr
            )?;


            // now that lateral connections are done, 
            // modify the heat transfer entity 

            self.pipe_fluid_array.set(fluid_array_clone.into())?;

            self.pipe_shell.set(pipe_shell_clone.into())?;

            self.insulation.set(insulation_array_clone.into())?;



        }
        // axial connections (insulation by default)
        // you can of course add new ones

        self.zero_power_bc_connection()?;

        Ok(())

    }


    /// for insulated pipes
    /// the end of each node should have a zero power boundary condition 
    /// connected to each of them at the bare minimum
    ///
    /// this function does exactly that
    ///
    /// to connect the rest of the heat transfer entities, 
    /// use the link to front or back methods within the 
    /// FluidArray or SolidColumn
    #[inline]
    fn zero_power_bc_connection(&mut self) -> Result<(),ThermalHydraulicsLibError>{

        let zero_power: Power = Power::ZERO;

        let mut zero_power_bc: HeatTransferEntity = 
        HeatTransferEntity::BoundaryConditions(
            BCType::UserSpecifiedHeatAddition(zero_power)
        );

        // constant heat addition interaction 

        let interaction: HeatTransferInteractionType = 
        HeatTransferInteractionType::UserSpecifiedHeatAddition;

        // now connect the twisted tape 


        self.pipe_fluid_array.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.pipe_fluid_array.link_to_back(&mut zero_power_bc,
            interaction)?;

        self.pipe_shell.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.pipe_shell.link_to_back(&mut zero_power_bc,
            interaction)?;

        self.insulation.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.insulation.link_to_back(&mut zero_power_bc,
            interaction)?;

        Ok(())
    }




    /// obtains air to insulation shell conductance
    ///
    /// it goes roughly to the middle of the insulation
    #[inline]
    pub fn get_ambient_surroundings_to_insulation_thermal_conductance(&mut self,
    h_air_to_pipe_surf: HeatTransfer) 
        -> Result<ThermalConductance,ThermalHydraulicsLibError> {
        // first, let's get a clone of the pipe_shell shell surface
        let mut insulation_clone: SolidColumn = 
        self.insulation.clone().try_into()?;

        let number_of_temperature_nodes = self.inner_nodes + 2;
        let heated_length = self.get_component_length();
        let insulation_id = self.insulation_id;
        let insulation_od = self.insulation_od;

        // next is to have pipe_shell inner conductance

        let insulation_shell_temperature: ThermodynamicTemperature 
        = insulation_clone.try_get_bulk_temperature()?;

        let cylinder_mid_diameter: Length = 0.5*(insulation_id+insulation_od);


        let node_length = heated_length / 
            number_of_temperature_nodes as f64;

        let pipe_air_conductance_interaction: HeatTransferInteractionType
        = HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(
                (insulation_clone.material_control_volume, 
                    (insulation_od-cylinder_mid_diameter).into(),
                    insulation_shell_temperature,
                    insulation_clone.pressure_control_volume),
                (h_air_to_pipe_surf,
                    insulation_od.into(),
                    node_length.into())
            );

        let pipe_air_nodal_thermal_conductance: ThermalConductance = try_get_thermal_conductance_based_on_interaction(
            self.ambient_temperature,
            insulation_shell_temperature,
            insulation_clone.pressure_control_volume,
            insulation_clone.pressure_control_volume,
            pipe_air_conductance_interaction,
        )?;


        return Ok(pipe_air_nodal_thermal_conductance);
    }


    /// obtains fluid_array to pipe_shell shell conductance
    #[inline]
    pub fn get_fluid_array_to_pipe_shell_conductance_no_wall_temp_correction(&mut self) 
        -> Result<ThermalConductance,ThermalHydraulicsLibError> {

        // the thermal conductance here should be based on the 
        // nusselt number correlation

        // before any calculations, I will first need a clone of 
        // the fluid_array fluid array and twisted tape array
        let mut fluid_array_clone: FluidArray = 
        self.pipe_fluid_array.clone().try_into()?;

        let mut pipe_shell_clone: SolidColumn = 
        self.pipe_shell.clone().try_into()?;

        // also need to get basic tmeperatures and mass flowrates 
        // only do this once because some of these methods involve 
        // cloning, which is computationally expensive

        let mass_flowrate: MassRate = 
        fluid_array_clone.get_mass_flowrate();

        let fluid_temperature: ThermodynamicTemperature 
        = fluid_array_clone.try_get_bulk_temperature()?;

        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let pipe_shell_surf_temperature: ThermodynamicTemperature 
        = pipe_shell_clone.try_get_bulk_temperature()?;

        let hydraulic_diameter = self.get_hydraulic_diameter();
        let flow_area: Area = self.get_cross_sectional_area_immutable();

        // firstly, reynolds 


        // flow area and hydraulic diameter are ok

        let mut fluid_array_clone: FluidArray = 
            self.pipe_fluid_array.clone().try_into()?;

        let fluid_material: LiquidMaterial
            = fluid_array_clone.material_control_volume.try_into()?;

        let solid_material: SolidMaterial
            = pipe_shell_clone.material_control_volume.try_into()?;

        let viscosity: DynamicViscosity = 
            fluid_material.try_get_dynamic_viscosity(fluid_temperature)?;

        // need to convert hydraulic diameter to an equivalent 
        // spherical diameter
        //
        // but for now, I'm going to use Re and Nu using hydraulic diameter 
        // and live with it for the time being
        //
        let reynolds_number: Ratio = 
            mass_flowrate/flow_area*hydraulic_diameter / viscosity;

        // next, bulk prandtl number 

        let bulk_prandtl_number: Ratio 
        = fluid_material.try_get_prandtl_liquid(
            fluid_temperature,
            atmospheric_pressure
        )?;


        //// surface prandtl number
        ////
        //let surface_prandtl_number: Ratio 
        //= LiquidMaterial::TherminolVP1.try_get_prandtl_liquid(
        //    pipe_shell_surf_temperature,
        //    atmospheric_pressure
        //).unwrap();
        ////note: we have an error here because fluid_array 
        // properties only range from 20 C to 180C,
        //
        // However, pipe_shell surface temperatures far exceed 180C 
        //
        // So the process will panic.
        // For now, we shall live within this temperature range

        // for this case, I will have the ciet heater nusselt 
        // number correlation
        //
        // constants are ignored, so we use the default method
        // and manually adjust the reynolds and prandtl numbers

        let mut pipe_prandtl_reynolds_data: GnielinskiData 
        = GnielinskiData::default();

        // no wall correction given for this case yet
        pipe_prandtl_reynolds_data.reynolds = reynolds_number;
        pipe_prandtl_reynolds_data.prandtl_bulk = bulk_prandtl_number;
        pipe_prandtl_reynolds_data.prandtl_wall = bulk_prandtl_number;
        pipe_prandtl_reynolds_data.length_to_diameter = 
            self.get_component_length_immutable()/
            self.get_hydraulic_diameter_immutable();

        // need to set darcy friction factor!
        // first let's get the fluid array out 


        let pipe_correlation = fluid_array_clone.fluid_component_loss_properties;

        // For this, the friction factor is the sum of pipe friction 
        // factor and form losses 
        // that's the convention
        let darcy_friction_factor: Ratio = pipe_correlation.darcy_friction_factor_fldk(
            reynolds_number)?;

        pipe_prandtl_reynolds_data.darcy_friction_factor = 
            darcy_friction_factor;

        let nusselt_estimate: Ratio = 
        pipe_prandtl_reynolds_data.get_nusselt_for_developing_flow()?;



        // now we can get the heat transfer coeff, 

        let h_to_fluid_array: HeatTransfer;

        let k_fluid_average: ThermalConductivity = 
        fluid_material.try_get_thermal_conductivity(
            fluid_temperature)?;

        h_to_fluid_array = nusselt_estimate * k_fluid_average / hydraulic_diameter;


        // and then get the convective resistance
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let heated_length = fluid_array_clone.get_component_length();
        let id = self.tube_id;
        let od = self.tube_od;


        let node_length = heated_length / 
            number_of_temperature_nodes as f64;


        // now I need to calculate resistance of the half length of the 
        // pipe_shell shell, which is an annular cylinder

        let cylinder_mid_diameter: Length = 0.5*(id+od);



        let fluid_array_pipe_shell_conductance_interaction: HeatTransferInteractionType
        = HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
                (solid_material.into(), 
                    (cylinder_mid_diameter - id).into(),
                    pipe_shell_surf_temperature,
                    atmospheric_pressure),
                (h_to_fluid_array,
                    id.into(),
                    node_length.into())
            );

        // now based on conductance interaction, 
        // we can obtain thermal conductance, the temperatures 
        // and pressures don't really matter
        //
        // this is because all the thermal conductance data 
        // has already been loaded into the thermal conductance 
        // interaction object

        let fluid_array_pipe_shell_nodal_thermal_conductance: ThermalConductance = 
            try_get_thermal_conductance_based_on_interaction(
                fluid_temperature,
                pipe_shell_surf_temperature,
                atmospheric_pressure,
                atmospheric_pressure,
                fluid_array_pipe_shell_conductance_interaction)?;


        return Ok(fluid_array_pipe_shell_nodal_thermal_conductance);
    }

    /// gets the reynolds number based on mass flworate and 
    /// temperature
    #[inline]
    pub fn get_reynolds_based_on_hydraulic_diameter_and_flow_area(
        &self,
        mass_flowrate: MassRate,
        temperature: ThermodynamicTemperature) -> Result<Ratio,ThermalHydraulicsLibError> {

        // flow area and hydraulic diameter are ok
        let flow_area: Area = self.get_cross_sectional_area_immutable();
        let hydraulic_diameter = self.get_hydraulic_diameter_immutable();

        let fluid_array_clone: FluidArray = 
            self.pipe_fluid_array.clone().try_into()?;

        let fluid_material: LiquidMaterial
            = fluid_array_clone.material_control_volume.try_into()?;

        let viscosity: DynamicViscosity = 
            fluid_material.try_get_dynamic_viscosity(temperature)?;

        // need to convert hydraulic diameter to an equivalent 
        // spherical diameter
        //
        // but for now, I'm going to use Re and Nu using hydraulic diameter 
        // and live with it for the time being
        //
        let reynolds: Ratio = 
        mass_flowrate/flow_area*hydraulic_diameter / viscosity;

        Ok(reynolds)

    }

    /// obtains fluid to pipe_shell conductance
    #[inline]
    pub fn get_fluid_node_pipe_shell_conductance(&mut self) 
        -> Result<ThermalConductance,ThermalHydraulicsLibError> {

        // the thermal conductance here should be based on the 
        // nusselt number correlation

        // before any calculations, I will first need a clone of 
        // the fluid array and pipe shell array
        let mut fluid_array_clone: FluidArray = 
        self.pipe_fluid_array.clone().try_into()?;

        let mut pipe_shell_clone: SolidColumn = 
        self.pipe_shell.clone().try_into().unwrap();

        // also need to get basic temperatures and mass flowrates 
        // only do this once because some of these methods involve 
        // cloning, which is computationally expensive

        let mass_flowrate: MassRate = 
        fluid_array_clone.get_mass_flowrate();

        let fluid_bulk_temperature: ThermodynamicTemperature 
        = fluid_array_clone.try_get_bulk_temperature()?;

        let solid_pressure = pipe_shell_clone.back_single_cv.pressure_control_volume;

        let pipe_shell_surf_temperature: ThermodynamicTemperature 
        = pipe_shell_clone.try_get_bulk_temperature()?;

        let hydraulic_diameter = 
        fluid_array_clone.get_hydraulic_diameter();

        let length = fluid_array_clone.get_component_length();

        // firstly, reynolds 

        let reynolds_number: Ratio = 
        self.get_reynolds_based_on_hydraulic_diameter_and_flow_area(
            mass_flowrate,
            fluid_bulk_temperature,
        )?;

        // materials
        
        let fluid_material: LiquidMaterial = 
            fluid_array_clone.back_single_cv.material_control_volume.try_into()?;

        let solid_material: SolidMaterial = 
            pipe_shell_clone.back_single_cv.material_control_volume.try_into()?;

        // next, bulk prandtl number 
        let bulk_prandtl_number: Ratio 
        = fluid_material.try_get_prandtl_liquid(
            fluid_bulk_temperature,
            solid_pressure
        )?;

        // surface prandtl number
        //
        let surface_prandtl_number: Ratio 
        = fluid_material.try_get_prandtl_liquid(
            pipe_shell_surf_temperature,
            solid_pressure
        )?;

        // for this case, I will have the Gnielinksi 
        // Correlation
        //
        // However, for that, I will need the length to diameter 
        // ratio, and the darcy_friction_factor


        let mut pipe_prandtl_reynolds_data: GnielinskiData 
        = GnielinskiData::default();

        pipe_prandtl_reynolds_data.reynolds = reynolds_number;
        pipe_prandtl_reynolds_data.prandtl_bulk = bulk_prandtl_number;
        pipe_prandtl_reynolds_data.prandtl_wall = surface_prandtl_number;

        pipe_prandtl_reynolds_data.darcy_friction_factor = 
            self.darcy_loss_correlation.darcy_friction_factor_fldk(
                reynolds_number).unwrap();

        pipe_prandtl_reynolds_data.length_to_diameter = 
        length/hydraulic_diameter;


        let heater_nusselt_correlation: NusseltCorrelation 
        =  NusseltCorrelation::PipeGnielinskiGeneric(
            pipe_prandtl_reynolds_data
        );

        let nusselt_estimate: Ratio = 
        heater_nusselt_correlation.try_get()?;

        // now we can get the heat transfer coeff, 

        let h: HeatTransfer;

        let k_fluid_average: ThermalConductivity = 
        fluid_material.try_get_thermal_conductivity(
            fluid_bulk_temperature)?;

        h = nusselt_estimate * k_fluid_average / hydraulic_diameter;

        // and then get the convective resistance
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let length = fluid_array_clone.get_component_length();
        let id = self.tube_id;
        let od = self.tube_od;

        let heat_transfer_area_total: Area = 
        length * id * PI;

        let heat_transfer_area_per_node: Area 
        = heat_transfer_area_total / 
        number_of_temperature_nodes as f64;

        let node_length = length / 
            number_of_temperature_nodes as f64;

        let fluid_to_pipe_shell_shell_average_conductance: ThermalConductance 
        = h * heat_transfer_area_per_node;

        let fluid_to_pipe_shell_shell_surface_node_resistance = 
        1.0/fluid_to_pipe_shell_shell_average_conductance;

        // now I need to calculate resistance of the half length of the 
        // pipe_shell shell, which is an annular cylinder

        let cylinder_mid_diameter: Length = 0.5*(id+od);

        let pipe_shell_conductivity = 
        solid_material.try_get_thermal_conductivity(
            pipe_shell_surf_temperature
        )?;

        let cylinder_node_conductance: ThermalConductance 
        = try_get_thermal_conductance_annular_cylinder(
            id,
            cylinder_mid_diameter,
            node_length,
            pipe_shell_conductivity
        )?;


        let cylinder_node_resistance = 
        1.0/cylinder_node_conductance;

        let cylinder_to_fluid_resistance = 
        cylinder_node_resistance + 
        fluid_to_pipe_shell_shell_surface_node_resistance;

        let cylinder_to_fluid_conductance: ThermalConductance 
        = 1.0/cylinder_to_fluid_resistance;

        return Ok(cylinder_to_fluid_conductance);
    }

    /// obtains pipe shell to insulation conductance
    #[inline]
    pub fn get_pipe_shell_to_insulation_conductance(
    &self) -> Result<ThermalConductance,ThermalHydraulicsLibError> {

        // first, make a clone of pipe shell and insulation

        let mut insulation_array_clone: SolidColumn = 
        self.insulation.clone().try_into()?;

        let mut pipe_shell_clone: SolidColumn = 
        self.pipe_shell.clone().try_into()?;


        // find the length of the array and node length

        let array_length =  self.get_component_length_immutable();

        let number_of_temperature_nodes = self.inner_nodes + 2;

        let node_length = array_length / 
        number_of_temperature_nodes as f64;

        // then we need to find the surface area of each node 
        // for steel to insulation_material, it will be 
        // the steel outer diameter or insulation inner_diameter
        
        let pipe_shell_mid_section_diameter = 0.5 * (self.tube_od 
        + self.tube_id);

        let insulation_material_mid_section_diameter = 0.5 * (self.insulation_id
        + self.insulation_od);

        let tube_od = self.tube_od;

        // next, thermal conductivities of both solid_pipe_material and insulation_material 

        let solid_pipe_material_shell_temperature = pipe_shell_clone.try_get_bulk_temperature() 
            ?;

        let solid_pipe_material: SolidMaterial = pipe_shell_clone.material_control_volume
            .try_into()?;

        let solid_pipe_material_conductivity: ThermalConductivity 
        = solid_pipe_material.try_get_thermal_conductivity(
            solid_pipe_material_shell_temperature
        )?;

        let insulation_material_shell_temperature = insulation_array_clone.try_get_bulk_temperature() 
            ?;

        let insulation_material: SolidMaterial = insulation_array_clone.material_control_volume
            .try_into()?;

        let insulation_material_conductivity: ThermalConductivity 
        = insulation_material.try_get_thermal_conductivity(
            insulation_material_shell_temperature
        )?;

        // we should be able to get the conductance now

        let insulation_material_layer_conductance: ThermalConductance = 
        try_get_thermal_conductance_annular_cylinder(
            tube_od,
            insulation_material_mid_section_diameter,
            node_length,
            insulation_material_conductivity
        )?;
        
        let solid_pipe_material_layer_conductance: ThermalConductance = 
        try_get_thermal_conductance_annular_cylinder(
            pipe_shell_mid_section_diameter,
            tube_od,
            node_length,
            solid_pipe_material_conductivity
        )?;

        // now that we have the conductances, we get the resistances 

        let insulation_material_resistance = 1.0/insulation_material_layer_conductance;
        let solid_pipe_material_resistance = 1.0/solid_pipe_material_layer_conductance;

        let total_resistance = insulation_material_resistance + solid_pipe_material_resistance;


        return Ok(1.0/total_resistance);
    }

    /// spawns a thread and moves the clone of the entire heater object into the 
    /// thread, "locking" it for parallel computation
    ///
    /// once that is done, the join handle is returned 
    /// which when unwrapped, returns the heater object
    pub fn lateral_connection_thread_spawn(&self,
    mass_flowrate: MassRate,
    heater_steady_state_power: Power) -> JoinHandle<Self>{

        let mut heater_clone = self.clone();

        // move ptr into a new thread 

        let join_handle = thread::spawn(
            move || -> Self {

                // carry out the connection calculations
                heater_clone.
                    lateral_and_miscellaneous_connections(
                        mass_flowrate,
                        heater_steady_state_power).unwrap();
                
                heater_clone

            }
        );

        return join_handle;

    }


}
