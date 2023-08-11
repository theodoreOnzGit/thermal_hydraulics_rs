use super::CartesianConduction1DArray;
use uom::si::f64::*;
use ndarray::*;
use uom::si::thermodynamic_temperature::kelvin;

impl CartesianConduction1DArray {

    /// gets bulk temperature of the array cv based on volume fraction 
    /// now, for solid and liquid, that would be sort of 
    /// a good approximation since the boussinesq approximation
    /// may work well for liquids
    ///
    #[inline]
    pub fn get_bulk_temperature(&mut self) -> 
    Result<ThermodynamicTemperature,String>{

        // for now, doing it quick and dirty, i'm going to obtain a volume 
        // averaged temperature 

        let volume_fraction_array_reference = 
        &self.volume_fraction_array;
        let temperature_array_reference = 
        &self.temperature_array_current_timestep;

        let mut vol_averaged_temperature_array_values: Array1<f64> 
        = Array::default(temperature_array_reference.len());

        for (idx, temperature_reference) in 
            temperature_array_reference.iter().enumerate() {
                //get the vol fraction 

                let vol_fraction: f64 = 
                volume_fraction_array_reference[idx];

                let vol_avg_temperature_component: f64
                = vol_fraction * (temperature_reference.get::<kelvin>());

                vol_averaged_temperature_array_values[idx] = 
                    vol_avg_temperature_component;

            }

        // sum it all up (these are float values) 

        let vol_averaged_temperature_kelvin: f64 
        = vol_averaged_temperature_array_values.sum();

        return Ok(ThermodynamicTemperature::new
            ::<kelvin>(vol_averaged_temperature_kelvin));


    }

    /// returns a clone of the temperature_array_current_timestep
    #[inline]
    pub fn get_temperature_vector(&mut self) -> 
    Result<Vec<ThermodynamicTemperature>, String> {
        let temp_array = self.temperature_array_current_timestep.clone();

        let mut temp_vector: Vec<ThermodynamicTemperature> = vec![];

        for temp_reference in temp_array.iter() {
            temp_vector.push(*temp_reference)
        }

        Ok(temp_vector)
    }
}
