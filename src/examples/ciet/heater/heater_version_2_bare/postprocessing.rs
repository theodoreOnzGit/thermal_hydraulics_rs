use uom::si::f64::*;

use super::HeaterVersion2Bare;

impl HeaterVersion2Bare {
    pub fn steel_shell_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.steel_shell.get_temperature_vector().unwrap()
    }

    pub fn _therminol_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.therminol_array.get_temperature_vector().unwrap()
    }

    pub fn _twisted_tape_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.twisted_tape_interior.get_temperature_vector().unwrap()
    }
}
