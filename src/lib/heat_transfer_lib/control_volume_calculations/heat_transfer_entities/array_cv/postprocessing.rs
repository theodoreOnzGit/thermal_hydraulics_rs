use uom::si::f64::*;
use crate::{heat_transfer_lib::{control_volume_calculations::heat_transfer_entities::ArrayCVType, thermophysical_properties::Material}, thermal_hydraulics_error::ThermalHydraulicsLibError};

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
                fluid_arr.try_get_bulk_temperature()
            },
            ArrayCVType::GenericColumn(solid_arr) => {
                solid_arr.try_get_bulk_temperature()
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
            ArrayCVType::GenericColumn(solid_arr) => {
                solid_arr.get_temperature_vector()
            },
        };

        temp_vector_result

        
    }

    /// gets pressure for array cv 
    pub fn get_array_cv_pressure(&mut self) -> Result<Pressure, ThermalHydraulicsLibError> {

        let pressure = match self {
            ArrayCVType::Cartesian1D(cv) => {
                cv.pressure_control_volume
            },
            ArrayCVType::GenericPipe(cv) => {
                cv.pressure_control_volume
            },
            ArrayCVType::GenericColumn(cv) => {
                cv.pressure_control_volume
            },
        };

        Ok(pressure)
    }

    /// gets pressure for array cv 
    pub fn get_array_cv_material(&mut self) -> Result<Material, ThermalHydraulicsLibError> {

        let material = match self {
            ArrayCVType::Cartesian1D(cv) => {
                cv.material_control_volume
            },
            ArrayCVType::GenericPipe(cv) => {
                cv.material_control_volume
            },
            ArrayCVType::GenericColumn(cv) => {
                cv.material_control_volume
            },
        };

        Ok(material)
    }

}
