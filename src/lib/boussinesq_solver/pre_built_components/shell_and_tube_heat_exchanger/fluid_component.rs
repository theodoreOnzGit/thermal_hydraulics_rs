use uom::si::f64::*;

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

    /// sets the tube side mass flowrate 
    pub fn set_tube_side_total_mass_flowrate(&mut self,
        mass_flowrate_over_all_tubes: MassRate) {

        let mut tube_side_fluid_array: FluidArray = 
        self.tube_side_parallel_fluid_array.clone().try_into().unwrap();


        let single_tube_mass_rate = 
            mass_flowrate_over_all_tubes / (self.number_of_tubes as f64);

        tube_side_fluid_array.set_mass_flowrate(single_tube_mass_rate);
        // unfortunately, this makes setting mass flowrate quite 
        // expensive as we need to clone it everytime

        self.tube_side_parallel_fluid_array.set(tube_side_fluid_array.into()).unwrap();

    }

    /// sets the tube side mass flowrate 
    pub fn set_shell_side_total_mass_flowrate(&mut self,
        mass_flowrate_through_shell: MassRate) {

        let mut shell_side_fluid_array: FluidArray = 
        self.shell_side_fluid_array.clone().try_into().unwrap();



        shell_side_fluid_array.set_mass_flowrate(mass_flowrate_through_shell);
        // unfortunately, this makes setting mass flowrate quite 
        // expensive as we need to clone it everytime

        self.shell_side_fluid_array.set(shell_side_fluid_array.into()).unwrap();

    }

}
