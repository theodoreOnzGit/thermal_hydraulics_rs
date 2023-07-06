use std::{error::Error, fmt};
/// generic thermophysical property error 
#[derive(Debug, Clone)]
pub struct ThermophysicalPropertyError;

impl Error for ThermophysicalPropertyError {}

impl fmt::Display for ThermophysicalPropertyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Error with calculating thermophysical_properties")
    }
}

/// this error happens if you use some of the private methods within 
/// the thermophysical_properties library 
#[derive(Debug, Clone)]
pub struct MaterialTypeError;

impl Error for MaterialTypeError {}

impl fmt::Display for MaterialTypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "You probably tried using a function of the wrong \n 
            phase type. For example, you want to calculate a liquid \n 
            enthalpy and inserted a solid material")
    }
}

/// this error happens if you use some of the private methods within 
/// the thermophysical_properties library 
#[derive(Debug, Clone)]
pub struct TemperatureRangeError;

impl Error for TemperatureRangeError { }

impl fmt::Display for TemperatureRangeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Temperature supplied for thermophysical_properties\n 
            functionw as out of range")
    }
}
