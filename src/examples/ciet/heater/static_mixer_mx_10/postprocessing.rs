use uom::si::f64::*;

use super::StaticMixerMX10;

impl StaticMixerMX10 {
    pub fn steel_shell_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.steel_shell.get_temperature_vector().unwrap()
    }

    pub fn therminol_array_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.therminol_array.get_temperature_vector().unwrap()
    }

    pub fn twisted_tape_temperature(&mut self) -> Vec<ThermodynamicTemperature>{
        self.insulation.get_temperature_vector().unwrap()
    }
}
