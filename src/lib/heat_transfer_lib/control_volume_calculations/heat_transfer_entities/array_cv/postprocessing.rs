use uom::si::f64::*;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;

/// contains implementations for ArrayCVType specific to post processing 
///
/// including getting bulk temperature
impl ArrayCVType {

    /// gets the bulk temperature for the ArrayCV
    pub fn get_bulk_temperature(&mut self) -> 
    Result<ThermodynamicTemperature,String>{

        match self {
            ArrayCVType::Cartesian1D(cartesian_1d_cv) => {
                cartesian_1d_cv.get_bulk_temperature()
            },
            ArrayCVType::GenericPipe(fluid_arr) => {
                Ok(fluid_arr.get_bulk_temperature().unwrap())
            },

        }

    }



}
