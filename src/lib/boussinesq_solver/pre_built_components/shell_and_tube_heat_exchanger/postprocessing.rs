use super::SimpleShellAndTubeHeatExchanger;
use uom::si::f64::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

impl SimpleShellAndTubeHeatExchanger {

    /// gets the temperature of the inner pipe shell array
    pub fn inner_pipe_shell_temperature(&mut self) -> 
        Result<Vec<ThermodynamicTemperature>, ThermalHydraulicsLibError>{
            self.inner_pipe_shell_array_for_single_tube.get_temperature_vector()
    }

    /// gets the temperature of the tube side fluid array
    pub fn inner_tube_fluid_array_temperature(&mut self) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.inner_pipe_shell_array_for_single_tube.get_temperature_vector()
    }

    /// gets the shell side fluid array temperature
    pub fn shell_side_fluid_array_temperature(&mut self,) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.shell_side_fluid_array.get_temperature_vector()
    }

    /// gets the shell side outer tube temperature 
    pub fn shell_side_outer_tube_array_temperature(&mut self,) -> 
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.outer_shell.get_temperature_vector()
    }

    /// gets the temperature of the insulation 
    pub fn insulation_array_temperature(&mut self,) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{

            if self.heat_exchanger_has_insulation {
                self.insulation_array.get_temperature_vector()
            } else {
                // if the heat exchanger is not even insulated, then 
                // there should be no insulation array
                todo!("simple Shell and Tube Heat exchanger is not insulated")
            }
    }

}
