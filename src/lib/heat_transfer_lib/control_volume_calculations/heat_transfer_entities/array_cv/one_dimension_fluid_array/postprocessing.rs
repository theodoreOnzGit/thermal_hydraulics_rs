use crate::fluid_mechanics_lib::prelude::FluidComponent;
use crate::heat_transfer_lib::thermophysical_properties::Material;
use crate::heat_transfer_lib::thermophysical_properties::prandtl::liquid_prandtl;
use crate::heat_transfer_lib::thermophysical_properties::
thermal_diffusivity::thermal_diffusivity;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use ndarray::*;
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::kelvin;
use super::FluidArray;
use uom::si::f64::*;

impl<const NUMBER_OF_NODES: usize> FluidArray<NUMBER_OF_NODES> {


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
