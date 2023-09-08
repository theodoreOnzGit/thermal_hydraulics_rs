use uom::si::f64::*;
use crate::{heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType, thermal_hydraulics_error::ThermalHydraulicsLibError};

/// contains implementations for ArrayCVType specific to post processing 
///
/// including getting bulk temperature
impl ArrayCVType {

    /// gets the bulk temperature for the ArrayCV
    pub fn get_bulk_temperature(&mut self) -> 
    Result<ThermodynamicTemperature,ThermalHydraulicsLibError>{

        match self {
            ArrayCVType::Cartesian1D(cartesian_1d_cv) => {
                cartesian_1d_cv.get_bulk_temperature()
            },
            ArrayCVType::GenericPipe(fluid_arr) => {
                Ok(fluid_arr.get_bulk_temperature().unwrap())
            },

        }

    }
    /// gets temperature vector for array cv
    pub fn get_temperature_vector(&mut self,) ->
    Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{

        let temp_vector_result = match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                match cartesian_array_cv.get_temperature_vector() {
                    Ok(vec) => Ok(vec),
                    Err(string) => Err(ThermalHydraulicsLibError::GenericStringError
                        (string)),
                }
            },
            ArrayCVType::GenericPipe(fluid_arr) => {
                fluid_arr.get_temperature_vector()
            },
        };

        temp_vector_result

        
    }



}
