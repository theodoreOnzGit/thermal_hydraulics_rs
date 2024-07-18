use std::f64::consts::PI;

use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::meter;

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

    /// returns tube side hydraulic diameter 
    /// by default, the internal diameter of the tube side
    #[inline]
    pub fn get_tube_side_hydraulic_diameter(&self) -> Length {
        return self.tube_side_id;
    }

    /// returns the shell side hydraulic diameter 
    /// this assumes that the shell side is a big tube with 
    /// a number of uniform circular tubes inside the big tube 
    ///
    /// from Du's paper, the formula here is:
    /// D_e = (D_i^2 - N_t d_o^2) / (D_i + N_t d_i)
    ///
    /// where 
    /// D_i  is shell_side_id
    /// d_i  is tube_side_id
    /// d_o is tube_side_od
    /// N_t is number_of_tubes
    #[inline]
    pub fn get_shell_side_hydraulic_diameter(&self) -> Length {


        // D_i 
        let shell_side_id = self.shell_side_id;
        
        //d_i 
        let tube_side_id = self.tube_side_id;

        //d_o 
        let tube_side_od = self.tube_side_od;

        // N_t 
        let number_of_tubes = self.number_of_tubes;

        // (D_i^2 - N_t d_o^2)

        let numerator: f64 = shell_side_id.get::<meter>().powf(2.0)
            - number_of_tubes as f64 
            * tube_side_od.get::<meter>().powf(2.0);
        
        // (D_i + N_t d_i)
        let denominator: f64 = 
            shell_side_id.get::<meter>() 
            + number_of_tubes as f64 * 
            tube_side_id.get::<meter>();

        // hydraulic diameter is in meters 

        let hydraulic_diameter: Length 
            = Length::new::<meter>(numerator/denominator);

        return hydraulic_diameter;
    }


    /// returns the shell side cross sectional area 
    /// this assumes that the shell side is a big tube with 
    /// a number of uniform circular tubes inside the big tube 
    ///
    /// from Du's paper, the formula here is:
    ///
    /// pi/4 * (D_i^2 - N_t d_o^2)
    pub fn get_shell_side_cross_sectional_area(&self) -> Area {

        // D_i 
        let shell_side_id = self.shell_side_id;
        
        //d_o 
        let tube_side_od = self.tube_side_od;

        // N_t 
        let number_of_tubes = self.number_of_tubes;

        // (D_i^2 - N_t d_o^2)

        let d_square_term: f64 = shell_side_id.get::<meter>().powf(2.0)
            - number_of_tubes as f64 
            * tube_side_od.get::<meter>().powf(2.0);
        
        let area_meter_sq_value: f64 
            = PI * 0.25 * d_square_term;

        let shell_side_xs_area: Area 
            = Area::new::<square_meter>(area_meter_sq_value);

        return shell_side_xs_area;
    }

}
