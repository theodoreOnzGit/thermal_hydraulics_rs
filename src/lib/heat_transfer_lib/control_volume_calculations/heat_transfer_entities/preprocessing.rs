use super::{HeatTransferEntity, CVType};
use uom::si::f64::*;

impl HeatTransferEntity {

    /// get maximum timestep 
    ///
    /// 
    pub fn get_max_timestep(
        entity: &mut HeatTransferEntity,
        max_temperature_change: TemperatureInterval) 
    -> Result<Time, String> {

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => 
                return Err("getting timestep not \n 
                    implemented for BoundaryConditions".to_string()),
        };

        let cv_timestep_result: 
        Result<Time, String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.get_max_timestep(max_temperature_change)
                },
                CVType::ArrayCV(cv) => {
                    cv.get_max_timestep(max_temperature_change)
                },
            };

        return cv_timestep_result;
    }
}
