use crate::boussinesq_solver::pre_built_components::heat_transfer_entities::HeatTransferEntity;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::{fluid_component_collection::fluid_component::FluidComponent, one_d_fluid_array_with_lateral_coupling::FluidArray};

use super::SimpleShellAndTubeHeatExchanger;

impl SimpleShellAndTubeHeatExchanger {


    /// clones the shell side fluid array, and converts it into a 
    /// fluid component
    pub fn get_clone_of_shell_side_fluid_component
        (&self) -> FluidComponent 
    {

        // first clone the heat transfer entity
        let shell_side_fluid_hte_clone: HeatTransferEntity 
            = self.shell_side_fluid_array.clone();

        // convert it into a fluid array
        let shell_side_fluid_array: FluidArray
            = shell_side_fluid_hte_clone.try_into().unwrap();

        return shell_side_fluid_array.into();

    }

}
