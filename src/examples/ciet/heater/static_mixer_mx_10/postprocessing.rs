use uom::si::f64::*;

use super::StaticMixerMX10;

impl StaticMixerMX10 {
    pub fn _steel_shell_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.steel_shell.get_temperature_vector().unwrap()
    }

    pub fn _therminol_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.therminol_array.get_temperature_vector().unwrap()
    }

    pub fn _insulation_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.insulation_array.get_temperature_vector().unwrap()
    }
}
