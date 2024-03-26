use uom::si::f64::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::SolidStructure;

impl SolidStructure {

    /// gets the temperature of the pipe shell array
    pub fn array_temperature(&mut self) -> 
        Result<Vec<ThermodynamicTemperature>, ThermalHydraulicsLibError>{
        self.solid_array.get_temperature_vector()
    }


}
