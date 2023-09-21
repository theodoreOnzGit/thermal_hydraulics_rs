use std::thread::{self, JoinHandle};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::{HeatTransferEntity, CVType};
use uom::si::f64::*;


impl HeatTransferEntity {

    /// for control volumes, this method allows you to 
    /// calculate the enthalpy of the next timestep and 
    /// set the 
    /// current timestep enthalpy as the enthalpy calculated  
    /// for the next timestep
    ///
    /// you are required to explicitly provide a timestep for this 
    pub fn advance_timestep(entity: &mut HeatTransferEntity,
    timestep: Time) -> Result<(), String> {

        // first match CV or BC, 
        // Boundary conditions don't need to advance timestep
        // so we can leave them be (it should return an Ok(()) value 
        // rather than an Err() value)

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => return Ok(()),
        };

        // once I have the cv_type enum, match it again

        let cv_advance_result: Result<(), String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.advance_timestep(timestep)
                },
                CVType::ArrayCV(cv) => {
                    cv.advance_timestep(timestep)
                },
            };

        return cv_advance_result;
    }


    /// for control volumes, this method allows you to 
    /// calculate the enthalpy of the next timestep and 
    /// set the 
    /// current timestep enthalpy as the enthalpy calculated  
    /// for the next timestep
    ///
    /// you are required to explicitly provide a timestep for this 
    pub fn advance_timestep_mut_self(&mut self,
    timestep: Time) -> Result<(), ThermalHydraulicsLibError> {

        // first match CV or BC, 
        // Boundary conditions don't need to advance timestep
        // so we can leave them be (it should return an Ok(()) value 
        // rather than an Err() value)

        let control_vol_type = match self {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => return Ok(()),
        };

        // once I have the cv_type enum, match it again

        let cv_advance_result: Result<(), ThermalHydraulicsLibError> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    Ok(single_cv.advance_timestep(timestep).unwrap())
                },
                CVType::ArrayCV(cv) => {
                    Ok(cv.advance_timestep(timestep).unwrap())
                },
            };

        return cv_advance_result;
    }

    /// spawns a handle to advance the timestep
    /// for parallel computation
    pub fn advance_timestep_mut_self_thread_spawn(&self,
        timestep: Time,) -> JoinHandle<Self> {

        // make a clone of the HeatTransferEntity (hte)
        let mut hte_clone = self.clone();

        // move ptr into a new thread 

        let join_handle = thread::spawn(
            move || -> Self {


                // carry out the connection calculations
                hte_clone.advance_timestep_mut_self(timestep).unwrap();
                
                hte_clone

            }
        );

        return join_handle;

    }

}
