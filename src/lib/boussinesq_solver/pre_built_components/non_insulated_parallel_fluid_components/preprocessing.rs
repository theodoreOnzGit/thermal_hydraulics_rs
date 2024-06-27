use std::thread::JoinHandle;
use std::thread;

use uom::ConstZero;
use uom::si::pressure::atmosphere;
use uom::si::f64::*;
use ndarray::*;
use super::NonInsulatedParallelFluidComponent;
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
impl NonInsulatedParallelFluidComponent {

    
    /// connects a heat transfer entity to the back of the parallel 
    /// fluid array given a heat transfer interaction
    pub fn connect_axially_to_back_of_parallel_fluid_array(
        &mut self,
        interaction: HeatTransferInteractionType,
        other_hte: HeatTransferEntity) -> Result<(), ThermalHydraulicsLibError>{

        // now, linking two heat transfer entities is fine, 
        // already coded that in. However, accounting for parallel tubes 
        // has not been done yet. This has to be hard coded
        //
        // unfortunately, hard coding it is also problematic because 
        // there are many many heat transfer interaction types
        //
        // so if there are n tubes, the power transferred to each tube (P) 
        // is now (P/n)
        //
        // this applies both for advection and conductance type interactions
        //
        // We could use the normal heat transfer interaction linking
        // for the other component, given the same total flowrate, the 
        // advection interaction does not really change 
        //
        // However, once the timestep advances, the heat transfer interaction 
        // for each single cv has to be a factor of 1/n
        //
        // so in fact, the simplest way to do parallel tubes is to alter the 
        // timestep advance step. At least for the axial parts
        //
        // Hard code it from here




        todo!()
    }

    /// connects a heat transfer entity to the front of the parallel 
    /// fluid array 
    pub fn connect_axially_to_front_of_parallel_fluid_array(
        &mut self,){

    }

    /// connects a boundary condition laterally to the shell 
    /// of the fluid array 
    pub fn connect_laterally_to_front_of_parallel_fluid_array(
        &mut self,){

    }

    /// lateral connections between shell and fluid arrays
    pub fn lateral_connection_shell_and_fluid_arrays(
        &mut self,){

    }

}
