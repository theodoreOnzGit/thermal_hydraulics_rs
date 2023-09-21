use uom::si::f64::*;

use super::StructuralSupport;

impl StructuralSupport {
    pub fn temperature_array(&mut self) -> Vec<ThermodynamicTemperature>{
        self.support_array.get_temperature_vector().unwrap()
    }

}
