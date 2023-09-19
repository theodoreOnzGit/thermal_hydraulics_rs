use super::StaticMixerMX10;
use thermal_hydraulics_rs::{prelude::alpha_nightly::*, heat_transfer_lib::{nusselt_correlations::{enums::NusseltCorrelation, input_structs::{NusseltPrandtlReynoldsData, GnielinskiData}}, control_volume_calculations::common_functions::try_get_thermal_conductance_annular_cylinder}};
use uom::{si::{area::square_inch, pressure::atmosphere}, ConstZero};
use ndarray::*;

impl StaticMixerMX10 {


    /// used to connect the arrays laterally 
    /// you'll need to set the mass flowrate and heater power
    ///
    /// executes serially, and uses lots of cloning, so it's 
    /// heavier in resource usage,
    ///
    /// unoptimised in this regard
    #[inline]
    pub fn start_lateral_connections(&mut self,
        mass_flowrate: MassRate,
        h_air_to_steel_surf: HeatTransfer){

        todo!();
        let heater_steady_state_power: Power = Power::ZERO;

        // first let's get all the conductances 

        let steel_to_air_conductance: ThermalConductance 
        = self.get_air_insulation_shell_conductance(
            h_air_to_steel_surf
        );

        self.set_mass_flowrate(mass_flowrate);

        let steel_surf_to_therminol_conductance: ThermalConductance 
        = self.get_therminol_node_steel_shell_conductance();

        let twisted_tape_to_therminol_conductance: ThermalConductance 
        = self.get_steel_to_fiberglass_conductance();

        // other stuff 
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let q_fraction_per_node: f64 = 1.0/ number_of_temperature_nodes as f64;
        let mut q_frac_arr: Array1<f64> = Array::default(number_of_temperature_nodes);
        q_frac_arr.fill(q_fraction_per_node);

        // then get the ambient temperature 

        let ambient_air_temp = self.ambient_temperature;

        // lateral connections 
        {
            // first i will need to create temperature vectors 

            let mut ambient_temperature_vector: Vec<ThermodynamicTemperature> 
            = Array1::default(number_of_temperature_nodes)
                .iter().map( |&temp| {
                    temp
                }
                ).collect();

            ambient_temperature_vector.fill(ambient_air_temp);

            let solid_temp_vector: Vec<ThermodynamicTemperature> 
            = self.steel_shell.get_temperature_vector().unwrap();

            let fluid_temp_vector: Vec<ThermodynamicTemperature> 
            = self.therminol_array.get_temperature_vector().unwrap();

            let twisted_tape_temp_vector: Vec<ThermodynamicTemperature> 
            = self.insulation_array.get_temperature_vector().unwrap();

            // clone each array and set them later

            let mut steel_shell_clone: SolidColumn = 
            self.steel_shell.clone().try_into() .unwrap();

            let mut therminol_array_clone: FluidArray = 
            self.therminol_array.clone().try_into().unwrap();

            let mut twisted_tape_array_clone: SolidColumn = 
            self.insulation_array.clone().try_into().unwrap();

            // second, fill them into the each array 
            
            // steel to air interaction

            steel_shell_clone.lateral_link_new_temperature_vector_avg_conductance(
                steel_to_air_conductance,
                ambient_temperature_vector
            ).unwrap();

            // steel shell to therminol interaction

            steel_shell_clone.lateral_link_new_temperature_vector_avg_conductance(
                steel_surf_to_therminol_conductance,
                fluid_temp_vector.clone()
            ).unwrap();

            therminol_array_clone.lateral_link_new_temperature_vector_avg_conductance(
                steel_surf_to_therminol_conductance,
                solid_temp_vector
            ).unwrap();

            // we also want to add a heat source to steel shell

            steel_shell_clone.lateral_link_new_power_vector(
                heater_steady_state_power,
                q_frac_arr
            ).unwrap();

            // now therminol to twisted tape interaction

            therminol_array_clone.
                lateral_link_new_temperature_vector_avg_conductance(
                    twisted_tape_to_therminol_conductance,
                    twisted_tape_temp_vector).unwrap();

            twisted_tape_array_clone. 
                lateral_link_new_temperature_vector_avg_conductance(
                    twisted_tape_to_therminol_conductance,
                    fluid_temp_vector).unwrap();

            // note, must set mass flowrate first 
            // otherwise there is by default zero flow through 
            // the array

            therminol_array_clone.set_mass_flowrate(
                mass_flowrate);


            // now that lateral connections are done, 
            // modify the heat transfer entity 

            self.therminol_array.set(therminol_array_clone.into()).unwrap();

            self.steel_shell.set(steel_shell_clone.into()).unwrap();

            self.insulation_array.set(twisted_tape_array_clone.into()
                ).unwrap();


        }
    }


    /// the end of each node should have a zero power boundary condition 
    /// connected to each of them at the bare minimum
    ///
    /// this function does exactly that
    ///
    /// to connect the rest of the heat transfer entities, 
    /// use the link to front or back methods within the 
    /// FluidArray or SolidColumn
    #[inline]
    fn zero_power_bc_connection(&mut self){

        let zero_power: Power = Power::ZERO;

        let mut zero_power_bc: HeatTransferEntity = 
        BCType::UserSpecifiedHeatAddition(zero_power).into();

        // constant heat addition interaction 

        let interaction: HeatTransferInteractionType = 
        HeatTransferInteractionType::UserSpecifiedHeatAddition;

        // now connect the twisted tape 

        self.insulation_array.link_to_back(&mut zero_power_bc,
            interaction).unwrap();


        self.insulation_array.link_to_front(&mut zero_power_bc,
            interaction).unwrap();

        self.therminol_array.link_to_front(&mut zero_power_bc,
            interaction).unwrap();

        self.therminol_array.link_to_back(&mut zero_power_bc,
            interaction).unwrap();

        self.steel_shell.link_to_front(&mut zero_power_bc,
            interaction).unwrap();

        self.steel_shell.link_to_back(&mut zero_power_bc,
            interaction).unwrap();
    }




    /// obtains air to steel shell conductance
    #[inline]
    pub fn get_air_insulation_shell_conductance(&mut self,
    h_air_to_insulation_surf: HeatTransfer) 
        -> ThermalConductance {

        // first, make a clone of steel 

        let mut fiberglass_clone: SolidColumn = 
        self.insulation_array.clone().try_into().unwrap();

        let mut therminol_clone: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        // find parameters for fiberglass conductance

        let number_of_temperature_nodes = self.inner_nodes + 2;
        let component_length: Length = therminol_clone.get_component_length();

        let fiberglass_shell_temperature = fiberglass_clone.try_get_bulk_temperature() 
            .unwrap();

        let fiberglass: SolidMaterial = fiberglass_clone.material_control_volume
            .try_into().unwrap();

        let fiberglass_conductivity: ThermalConductivity 
        = fiberglass.try_get_thermal_conductivity(
            fiberglass_shell_temperature
        ).unwrap();

        let node_length: Length  = component_length / 
        number_of_temperature_nodes as f64;

        let od: Length = self.insulation_outer_diameter;

        let insulation_mid_length: Length 
        = 0.5 * (od + self.insulation_inner_diameter);

        let fiberglass_layer_conductance: ThermalConductance = 
        try_get_thermal_conductance_annular_cylinder(
            insulation_mid_length,
            od,
            node_length,
            fiberglass_conductivity
        ).unwrap();

        // find parameters for conductance due to convection 
        // resistance 

        let node_area: Area = od * PI * node_length;

        let air_convection_conductance: ThermalConductance
        = node_area * h_air_to_insulation_surf;

        let total_resistance = 
        1.0/air_convection_conductance + 
        1.0/fiberglass_layer_conductance;

        return 1.0/total_resistance;

    }


    /// obtains therminol to steel shell conductance
    #[inline]
    pub fn get_therminol_node_steel_shell_conductance(&mut self) 
        -> ThermalConductance {

        // the thermal conductance here should be based on the 
        // nusselt number correlation

        // before any calculations, I will first need a clone of 
        // the therminol fluid array and twisted tape array
        let mut therminol_fluid_array_clone: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        let mut steel_shell_clone: SolidColumn = 
        self.steel_shell.clone().try_into().unwrap();

        // also need to get basic tmeperatures and mass flowrates 
        // only do this once because some of these methods involve 
        // cloning, which is computationally expensive

        let mass_flowrate: MassRate = 
        therminol_fluid_array_clone.get_mass_flowrate();

        let bulk_temperature: ThermodynamicTemperature 
        = therminol_fluid_array_clone.try_get_bulk_temperature().unwrap();

        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let steel_surf_temperature: ThermodynamicTemperature 
        = steel_shell_clone.try_get_bulk_temperature().unwrap();

        let hydraulic_diameter = 
        therminol_fluid_array_clone.get_hydraulic_diameter();

        let length = therminol_fluid_array_clone.get_component_length();

        // firstly, reynolds 

        let reynolds_number: Ratio = 
        self.mx10_hydraulic_diameter_reynolds(
            mass_flowrate,
            bulk_temperature,
        );

        // next, bulk prandtl number 

        let bulk_prandtl_number: Ratio 
        = LiquidMaterial::TherminolVP1.try_get_prandtl_liquid(
            bulk_temperature,
            atmospheric_pressure
        ).unwrap();

        // surface prandtl number
        //
        let surface_prandtl_number: Ratio 
        = LiquidMaterial::TherminolVP1.try_get_prandtl_liquid(
            steel_surf_temperature,
            atmospheric_pressure
        ).unwrap();

        // for this case, I will have the Gnielinksi 
        // Correlation
        //
        // However, for that, I will need the length to diameter 
        // ratio, and the darcy_friction_factor


        let mut mx10_prandtl_reynolds_data: GnielinskiData 
        = GnielinskiData::default();

        mx10_prandtl_reynolds_data.reynolds = reynolds_number;
        mx10_prandtl_reynolds_data.prandtl_bulk = bulk_prandtl_number;
        mx10_prandtl_reynolds_data.prandtl_wall = surface_prandtl_number;

        mx10_prandtl_reynolds_data.darcy_friction_factor = 
            self.darcy_loss_correlation.darcy_friction_factor_fldk(
                reynolds_number).unwrap();

        mx10_prandtl_reynolds_data.length_to_diameter = 
        length/hydraulic_diameter;


        let heater_nusselt_correlation: NusseltCorrelation 
        =  NusseltCorrelation::PipeGnielinskiGeneric(
            mx10_prandtl_reynolds_data
        );

        let nusselt_estimate: Ratio = 
        heater_nusselt_correlation.try_get().unwrap();

        // now we can get the heat transfer coeff, 

        let h: HeatTransfer;

        let k_fluid_average: ThermalConductivity = 
        LiquidMaterial::TherminolVP1.try_get_thermal_conductivity(
            bulk_temperature).unwrap();

        h = nusselt_estimate * k_fluid_average / hydraulic_diameter;

        // and then get the convective resistance
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let length = therminol_fluid_array_clone.get_component_length();
        let id = Length::new::<meter>(0.0381);
        let od = Length::new::<meter>(0.04);

        let heat_transfer_area_total: Area = 
        length * id * PI;

        let heat_transfer_area_per_node: Area 
        = heat_transfer_area_total / 
        number_of_temperature_nodes as f64;

        let node_length = length / 
            number_of_temperature_nodes as f64;

        let therminol_to_steel_shell_average_conductance: ThermalConductance 
        = h * heat_transfer_area_per_node;

        let therminol_to_steel_shell_surface_node_resistance = 
        1.0/therminol_to_steel_shell_average_conductance;

        // now I need to calculate resistance of the half length of the 
        // steel shell, which is an annular cylinder

        let cylinder_mid_diameter: Length = 0.5*(id+od);

        let steel_conductivity: ThermalConductivity = 
        SolidMaterial::SteelSS304L.try_get_thermal_conductivity(
            steel_surf_temperature
        ).unwrap();

        let cylinder_node_conductance: ThermalConductance 
        = try_get_thermal_conductance_annular_cylinder(
            id,
            cylinder_mid_diameter,
            node_length,
            steel_conductivity
        ).unwrap();


        let cylinder_node_resistance = 
        1.0/cylinder_node_conductance;

        let cylinder_to_therminol_resistance = 
        cylinder_node_resistance + 
        therminol_to_steel_shell_surface_node_resistance;

        let cylinder_to_therminol_conductance: ThermalConductance 
        = 1.0/cylinder_to_therminol_resistance;

        return cylinder_to_therminol_conductance;
    }

    #[inline]
    pub fn mx10_hydraulic_diameter_reynolds(&mut self, 
        mass_flowrate: MassRate,
        temperature: ThermodynamicTemperature) -> Ratio {

        let flow_area: Area = self.flow_area;
        let hydraulic_diameter = self.get_hydraulic_diameter();
        let viscosity: DynamicViscosity = 
        LiquidMaterial::TherminolVP1.try_get_dynamic_viscosity(
            temperature).unwrap();

        // need to convert hydraulic diameter to an equivalent 
        // spherical diameter
        //
        // but for now, I'm going to use Re and Nu using hydraulic diameter 
        // and live with it for the time being
        //
        let reynolds: Ratio = 
        mass_flowrate/flow_area*hydraulic_diameter / viscosity;

        reynolds

    }


    /// obtains therminol to twisted tape conductance 
    /// based on approx wakao correlation
    #[inline]
    pub fn get_steel_to_fiberglass_conductance(
    &self) -> ThermalConductance {

        // first, make a clone of steel and fiberglass nodes 

        let mut fiberglass_clone: SolidColumn = 
        self.insulation_array.clone().try_into().unwrap();

        let mut steel_clone: SolidColumn = 
        self.steel_shell.clone().try_into().unwrap();

        let mut therminol_clone: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        // find the length of the array and node length

        let array_length =  therminol_clone.get_component_length();

        let number_of_temperature_nodes = self.inner_nodes + 2;

        let node_length = array_length / 
        number_of_temperature_nodes as f64;

        // then we need to find the surface area of each node 
        // for steel to fiberglass, it will be 
        // the steel outer diameter or insulation inner_diameter
        
        let steel_mid_section_diameter = 0.5 * (self.tube_outer_diameter 
        + self.tube_inner_diameter);

        let fiberglass_mid_section_diameter = 0.5 * (self.insulation_inner_diameter
        + self.insulation_outer_diameter);

        let steel_od = self.tube_outer_diameter;

        // next, thermal conductivities of both steel and fiberglass 

        let steel_shell_temperature = steel_clone.try_get_bulk_temperature() 
            .unwrap();

        let steel: SolidMaterial = steel_clone.material_control_volume
            .try_into().unwrap();

        let steel_conductivity: ThermalConductivity 
        = steel.try_get_thermal_conductivity(
            steel_shell_temperature
        ).unwrap();

        let fiberglass_shell_temperature = fiberglass_clone.try_get_bulk_temperature() 
            .unwrap();

        let fiberglass: SolidMaterial = fiberglass_clone.material_control_volume
            .try_into().unwrap();

        let fiberglass_conductivity: ThermalConductivity 
        = fiberglass.try_get_thermal_conductivity(
            fiberglass_shell_temperature
        ).unwrap();

        // we should be able to get the conductance now

        let fiberglass_layer_conductance: ThermalConductance = 
        try_get_thermal_conductance_annular_cylinder(
            steel_od,
            fiberglass_mid_section_diameter,
            node_length,
            fiberglass_conductivity
        ).unwrap();
        
        let steel_layer_conductance: ThermalConductance = 
        try_get_thermal_conductance_annular_cylinder(
            steel_mid_section_diameter,
            steel_od,
            node_length,
            steel_conductivity
        ).unwrap();

        // now that we have the conductances, we get the resistances 

        let fiberglass_resistance = 1.0/fiberglass_layer_conductance;
        let steel_resistance = 1.0/steel_layer_conductance;

        let total_resistance = fiberglass_resistance + steel_resistance;


        return 1.0/total_resistance;
    }
}


