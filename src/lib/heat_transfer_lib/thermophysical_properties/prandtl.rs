use uom::si::f64::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use super::Material;
use super::dynamic_viscosity::dynamic_viscosity;
use super::specific_heat_capacity::specific_heat_capacity;
use super::thermal_conductivity::thermal_conductivity;


/// provides the prandtl number
#[inline]
pub fn liquid_prandtl(material: Material,
    temperature: ThermodynamicTemperature,
    pressure: Pressure) -> Result<Ratio, ThermalHydraulicsLibError> {

    // get mu 
    let mu: DynamicViscosity = dynamic_viscosity(
        material,
        temperature,
        pressure)?;

    // get cp 
    let cp: SpecificHeatCapacity = specific_heat_capacity(
        material,
        temperature,
        pressure)?;

    // get k 
    let k: ThermalConductivity = thermal_conductivity(
        material,
        temperature,
        pressure)?;

    Ok(mu*cp/k)

}
