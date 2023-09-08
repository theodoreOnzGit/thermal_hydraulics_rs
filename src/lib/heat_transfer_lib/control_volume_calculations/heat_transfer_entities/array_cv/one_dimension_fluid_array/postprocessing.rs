use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::FluidArray;
use uom::si::f64::*;

impl FluidArray {


    /// obtains a clone of the temperature vector within the CV 
    /// thus obtaining the temperature profile
    pub fn get_temperature_vector(&self) -> Result<
    Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        let mut temperature_vec: Vec<ThermodynamicTemperature> = vec![];

        for temperature in self.temperature_array_current_timestep.iter() {
            temperature_vec.push(*temperature);
        }

        return Ok(temperature_vec);
    }

}
