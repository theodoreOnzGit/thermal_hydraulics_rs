use uom::si::f64::*;

use super::StructuralSupport;

impl StructuralSupport {
    pub fn steel_shell_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.steel_shell.get_temperature_vector().unwrap()
    }

    pub fn therminol_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.therminol_array.get_temperature_vector().unwrap()
    }

    pub fn insulation_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.insulation_array.get_temperature_vector().unwrap()
    }
}
