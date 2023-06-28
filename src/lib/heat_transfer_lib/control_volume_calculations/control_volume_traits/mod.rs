//! in this module, I want to house traits for control volumes
//!
//! a control volume in the heat transfer context will just
//! be a fixed region of space.
//!
//! This region will have fluid flowing in and out of the control volume
//! thus bringing enthalpy in and out of the control volume
//!
//! also there is heat and work just like the
//! first law of thermodynamics
//!
//! the idea here is to calculate the enthalpy of
//! the control volume at the next timestep so to speak
//!
//! the thing here is that the control volume is perfectly well
//! mixed so that the outlet temperature is equal to the bulk temperature
//! of the control volume
//!

use uom::si::f64::*;
use super::common_functions::*;

/// traits for pipe like control volumes
/// ie, one flow in, one flow out
mod pipe;
//use pipe::*;


/// this is the primary trait for control volume,
/// meant for making trait objects
trait ControlVolume: ControlVolumeAssociatedFunctions {

    // for each calculation, we get our control volume enthalpy

    /// gets the control volume specific enthalpy
    fn get_control_volume_specific_enthalpy(&self) -> AvailableEnergy;

    /// gets the control volume mass
    fn get_control_volume_mass(&self) -> Mass;

    /// get sum of enthalpy flows in
    fn get_sum_of_enthalpy_flows_in(&self) -> Power;

    /// get sum of enthalpy flows out
    fn get_sum_of_enthalpy_flows_out(&self) -> Power;

    /// get_heat_supplied_to_system
    fn get_heat_supplied_to_system(&self) -> Power;

    /// get work done on system
    fn get_work_done_on_system(&self) -> Power;

    /// convert enthalpy to temperature
    fn get_fluid_temperature_from_specific_enthalpy(
        &self,
        fluid_enthalpy: AvailableEnergy) -> ThermodynamicTemperature;


    /// calculate next timestep control volume temperature
    ///
    /// I leave timestep here as an input because
    /// we may want to have adjustable timestep
    fn calculate_control_vol_temp_next_timestep(
        &self,
        timestep: Time) -> ThermodynamicTemperature{

        let current_timestep_control_volume_specific_enthalpy 
            = self.get_control_volume_specific_enthalpy();

        let current_timestep_control_volume_mass  
            = self.get_control_volume_mass();

        let enthalpy_in 
            = self.get_sum_of_enthalpy_flows_in();

        let enthalpy_out 
            = self.get_sum_of_enthalpy_flows_out();

        let heat_supplied_to_system 
            = self.get_heat_supplied_to_system();

        let work_done_on_system
            = self.get_work_done_on_system();

        // we may want to change this in future to involve Courant
        // number calcs
        let next_timestep_specific_enthalpy
            = Self::calculate_specific_enthalpy_at_next_timestep(
                current_timestep_control_volume_specific_enthalpy, 
                current_timestep_control_volume_mass, 
                timestep, 
                enthalpy_out, 
                enthalpy_in, 
                heat_supplied_to_system, 
                work_done_on_system);



        let next_timestep_control_vol_temperature
            = self.get_fluid_temperature_from_specific_enthalpy(
                next_timestep_specific_enthalpy.unwrap());

        return next_timestep_control_vol_temperature;

    }

    
}

/// contains associated functions for the control
/// volume,
///
/// not meant for making into trait objects
trait ControlVolumeAssociatedFunctions {

    /// calculates enthalpy of a control volume at the next timestep
    ///
    /// this assumes that the mass of the control volume does not change
    /// in any significant way
    ///
    /// this function does not check for Courant number
    fn calculate_specific_enthalpy_at_next_timestep(
        current_timestep_control_volume_specific_enthalpy: AvailableEnergy,
        current_timestep_control_volume_mass: Mass,
        timestep: Time,
        enthalpy_out: Power,
        enthalpy_in: Power,
        heat_supplied_to_system: Power,
        work_done_on_system: Power,
        ) -> Result<AvailableEnergy,AvailableEnergy> {

        let control_volume_enthalpy_current_timestep: Energy =
            current_timestep_control_volume_mass
            * current_timestep_control_volume_specific_enthalpy;

        let control_volume_enthalpy_next_timestep = 
            get_control_volume_enthalpy_next_timestep(
                timestep, 
                enthalpy_out, 
                enthalpy_in, 
                heat_supplied_to_system, 
                work_done_on_system, 
                control_volume_enthalpy_current_timestep);

        // we are assuming control volume mass does not change
        let specific_enthalpy_next_timestep = 
            control_volume_enthalpy_next_timestep/
            current_timestep_control_volume_mass;

        return Ok(specific_enthalpy_next_timestep);
    }
    

}
