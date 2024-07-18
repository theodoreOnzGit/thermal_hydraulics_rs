use std::thread::JoinHandle;
use std::thread;

use uom::ConstZero;
use uom::si::pressure::atmosphere;
use uom::si::f64::*;
use ndarray::*;
use super::SimpleShellAndTubeHeatExchanger;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component::FluidComponent;
use crate::boussinesq_solver::{heat_transfer_correlations::nusselt_number_correlations::input_structs::GnielinskiData, pre_built_components::heat_transfer_entities::preprocessing::try_get_thermal_conductance_based_on_interaction};
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boundary_conditions::BCType;
use crate::boussinesq_solver::pre_built_components::heat_transfer_entities::HeatTransferEntity;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidComponentTrait;

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

// preprocessing is where heat transfer entities 
// are connected to each other whether axially or laterally
//
// it is important to pay attention to the code here because it contains 
// logic of how heat transfer entites interact with the arrayCVs within 
// this parallel fluid component 
//
// first, we want a function to connect the components via heat transfer 
// interactions axially
//
// then we need to connect the parallel fluid arrays laterally to some 
// boundary condition
impl SimpleShellAndTubeHeatExchanger {

    #[inline]
    pub fn lateral_and_miscellaneous_connections(&mut self,
    ) -> Result<(), ThermalHydraulicsLibError>
    {

        // first let's get all the conductances 
        let heat_transfer_to_ambient = self.heat_transfer_to_ambient;

        let outer_node_to_air_conductance = 
            self.get_air_to_single_shell_nodal_shell_conductance(
                heat_transfer_to_ambient)?;

        let insulation_to_shell_conductance: ThermalConductance;

        let shell_to_shell_side_fluid_conductance: ThermalConductance;

        // for the parallel tube bundle, we have to be extra careful 
        // because we are accounting for multiple tubes 
        //
        // for the tube, it is one single pipe so to speak, 
        // Therefore, the tube fluid to tube metallic side should be 
        // for a single tube only 
        //
        // However, for shell side fluid to tube, we could use the 
        // conductance to a single tube or to the tube bundle
        // For clarity, I just named the conductance by the tube bundle.
        //
        // However, I also define a single tube to shell side fluid conductance
        // This avoids ambiguity when dealing with the conductance arrays
        //
        let tube_bundle_to_shell_side_fluid_conductance: ThermalConductance;
        let single_tube_to_shell_side_fluid_conductance: ThermalConductance;
        let single_tube_to_tube_side_fluid_conductance: ThermalConductance;

        // axial connections  (adiabatic by default)
        self.zero_power_bc_axial_connection()?;

        todo!()

    }


    /// the end of each node should have a zero power boundary condition 
    /// connected to each of them at the bare minimum
    ///
    /// this function does exactly that
    ///
    /// to connect the rest of the heat transfer entities, 
    /// use the link to front or back methods within the 
    /// FluidArrays or SolidColumns
    #[inline]
    fn zero_power_bc_axial_connection(&mut self) -> Result<(),ThermalHydraulicsLibError>{

        let zero_power: Power = Power::ZERO;

        let mut zero_power_bc: HeatTransferEntity = 
        HeatTransferEntity::BoundaryConditions(
            BCType::UserSpecifiedHeatAddition(zero_power)
        );

        // constant heat addition interaction 

        let interaction: HeatTransferInteractionType = 
        HeatTransferInteractionType::UserSpecifiedHeatAddition;

        // now connect the four or five arrays 
        // two solid arrays (or three if insulation is switched on) 
        // and two fluid arrays


        // tube side
        self.tube_side_parallel_fluid_array.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.tube_side_parallel_fluid_array.link_to_back(&mut zero_power_bc,
            interaction)?;

        self.pipe_shell.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.pipe_shell.link_to_back(&mut zero_power_bc,
            interaction)?;

        // shell side
        self.shell_side_fluid_array.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.shell_side_fluid_array.link_to_back(&mut zero_power_bc,
            interaction)?;

        self.outer_shell.link_to_front(&mut zero_power_bc,
            interaction)?;

        self.outer_shell.link_to_back(&mut zero_power_bc,
            interaction)?;

        // insulation 

        if self.heat_exchanger_has_insulation {

            self.insulation_array.link_to_front(&mut zero_power_bc,
                interaction)?;

            self.insulation_array.link_to_back(&mut zero_power_bc,
                interaction)?;

        }


        Ok(())
    }


    /// obtains air to shell and tube heat exchanger 
    /// outer array conductance 
    ///
    /// The outer array will be insulation if insulation is switched on,
    /// or the outer shell if insulation is switched off
    #[inline]
    pub fn get_air_to_single_shell_nodal_shell_conductance(&mut self,
        h_air_to_pipe_surf: HeatTransfer) 
        -> Result<ThermalConductance,ThermalHydraulicsLibError> 
    {

        // for conductance calculations (no radiation), 
        // it is important to get the temperatures of the ambient 
        // surroundings and the dimensions of the outer shell or insulation

        let heated_length: Length;
        let id: Length;
        let od: Length;
        let outer_node_temperature: ThermodynamicTemperature;
        // shell and tube heat excanger (STHE) to air interaction
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let outer_solid_array_clone: SolidColumn;

        if self.heat_exchanger_has_insulation {
            // if there's insulation, the id is the outer diameter of 
            // the shell. 

            id = self.shell_side_od;
            od = self.shell_side_od + self.insulation_thickness;

            // heated length is the shell side length 
            // first I need the fluid array as a fluid component

            let shell_side_fluid_component_clone: FluidComponent 
                = self.get_clone_of_shell_side_fluid_component();

            // then i need to get the component length 
            heated_length = shell_side_fluid_component_clone
                .get_component_length_immutable();

            // surface temperature is the insulation bulk temperature 
            // (estimated)

            let mut shell_side_fluid_array: FluidArray = 
                shell_side_fluid_component_clone.try_into().unwrap();

            outer_node_temperature = shell_side_fluid_array
                .try_get_bulk_temperature()?;

            // the outer node clone is insulation if it is switched on
            outer_solid_array_clone = 
                self.insulation_array.clone().try_into()?;

        } else {
            // if there's no insulation, the id is the inner diameter of 
            // the shell
            // od is outer diameter of the shell

            id = self.shell_side_id;
            od = self.shell_side_od;

            // heated length is the shell side length 
            // first I need the fluid array as a fluid component

            let shell_side_fluid_component_clone: FluidComponent 
                = self.get_clone_of_shell_side_fluid_component();

            // then i need to get the component length 
            heated_length = shell_side_fluid_component_clone
                .get_component_length_immutable();

            // surface temperature is the insulation bulk temperature 
            // (estimated)

            let mut shell_side_fluid_array: FluidArray = 
                shell_side_fluid_component_clone.try_into().unwrap();

            outer_node_temperature = shell_side_fluid_array
                .try_get_bulk_temperature()?;

            // the outer node clone is outer shell array if it is switched off
            outer_solid_array_clone = 
                self.outer_shell.clone().try_into()?;

        }

        let cylinder_mid_diameter: Length = 0.5*(id+od);


        let node_length = heated_length / 
            number_of_temperature_nodes as f64;

        let outer_node_air_conductance_interaction: HeatTransferInteractionType
        = HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(
                (outer_solid_array_clone.material_control_volume, 
                    (od-cylinder_mid_diameter).into(),
                    outer_node_temperature,
                    outer_solid_array_clone.pressure_control_volume),
                (h_air_to_pipe_surf,
                    od.into(),
                    node_length.into())
            );

        let outer_node_air_nodal_thermal_conductance: ThermalConductance = try_get_thermal_conductance_based_on_interaction(
            self.ambient_temperature,
            outer_node_temperature,
            outer_solid_array_clone.pressure_control_volume,
            outer_solid_array_clone.pressure_control_volume,
            outer_node_air_conductance_interaction,
        )?;


        return Ok(outer_node_air_nodal_thermal_conductance);
    }


    /// obtains fluid to pipe shell conductance
    #[inline]
    pub fn get_single_tube_side_fluid_array_node_to_pipe_shell_conductance_with_wall_temp_correction(
        &mut self) 
        -> Result<ThermalConductance,ThermalHydraulicsLibError> 
    {

        // the thermal conductance here should be based on the 
        // nusselt number correlation

        // before any calculations, I will first need a clone of 
        // the fluid array and twisted tape array
        let mut fluid_array_clone: FluidArray = 
            self.tube_side_parallel_fluid_array.clone().try_into()?;

        let mut pipe_shell_clone: SolidColumn = 
            self.pipe_shell.clone().try_into()?;

        // also need to get basic temperatures and mass flowrates 
        // only do this once because some of these methods involve 
        // cloning, which is computationally expensive

        let single_tube_mass_flowrate: MassRate = 
            fluid_array_clone.get_mass_flowrate();

        let fluid_temperature: ThermodynamicTemperature 
            = fluid_array_clone.try_get_bulk_temperature()?;

        let wall_temperature: ThermodynamicTemperature 
            = pipe_shell_clone.try_get_bulk_temperature()?;

        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let pipe_shell_surf_temperature: ThermodynamicTemperature 
            = pipe_shell_clone.try_get_bulk_temperature()?;

        let single_tube_hydraulic_diameter = 
            fluid_array_clone.get_hydraulic_diameter();
        let bundled_tubes_flow_area: Area = 
            fluid_array_clone.get_cross_sectional_area_immutable();

        let single_tube_flow_area: Area = 
            bundled_tubes_flow_area / (self.number_of_tubes as f64);


        // flow area and hydraulic diameter are ok


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
        let reynolds_number_single_tube: Ratio = 
            single_tube_mass_flowrate/
            single_tube_flow_area
            *single_tube_hydraulic_diameter / viscosity;

        // next, bulk prandtl number 

        let bulk_prandtl_number: Ratio 
            = fluid_material.try_get_prandtl_liquid(
                fluid_temperature,
                atmospheric_pressure
            )?;

        // then wall prandtl number

        let wall_prandtl_number: Ratio 
            = fluid_material.try_get_prandtl_liquid(
                wall_temperature,
                atmospheric_pressure
            )?;


        let mut pipe_prandtl_reynolds_data: GnielinskiData 
            = GnielinskiData::default();

        // wall correction given for this case yet
        pipe_prandtl_reynolds_data.reynolds = reynolds_number_single_tube;
        pipe_prandtl_reynolds_data.prandtl_bulk = bulk_prandtl_number;
        pipe_prandtl_reynolds_data.prandtl_wall = wall_prandtl_number;
        pipe_prandtl_reynolds_data.length_to_diameter = 
            fluid_array_clone.get_component_length_immutable()/
            fluid_array_clone.get_hydraulic_diameter_immutable();

        // I need to use Nusselt correlations present in this struct 
        //
        // wall correction is done here
        //
        // this uses the gnielinski correlation for pipes or tubes

        let nusselt_estimate = 
            pipe_prandtl_reynolds_data.
            get_nusselt_for_developing_flow()?;



        // now we can get the heat transfer coeff, 

        let h_to_fluid: HeatTransfer;

        let k_fluid_average: ThermalConductivity = 
            fluid_material.try_get_thermal_conductivity(
                fluid_temperature)?;

        h_to_fluid = nusselt_estimate * k_fluid_average / single_tube_hydraulic_diameter;


        // and then get the convective resistance
        let number_of_temperature_nodes = self.inner_nodes + 2;
        let heated_length = fluid_array_clone.get_component_length();
        let id = self.tube_side_id;
        let od = self.tube_side_od;


        let node_length = heated_length / 
            number_of_temperature_nodes as f64;


        // now I need to calculate resistance of the half length of the 
        // pipe shell, which is an annular cylinder

        let cylinder_mid_diameter: Length = 0.5*(id+od);



        let fluid_pipe_shell_conductance_interaction: HeatTransferInteractionType
            = HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
                (solid_material.into(), 
                 (cylinder_mid_diameter - id).into(),
                 pipe_shell_surf_temperature,
                 atmospheric_pressure),
                 (h_to_fluid,
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

        let fluid_pipe_shell_nodal_thermal_conductance: ThermalConductance = 
            try_get_thermal_conductance_based_on_interaction(
                fluid_temperature,
                pipe_shell_surf_temperature,
                atmospheric_pressure,
                atmospheric_pressure,
                fluid_pipe_shell_conductance_interaction)?;


        return Ok(fluid_pipe_shell_nodal_thermal_conductance);
    }
}
