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



}
