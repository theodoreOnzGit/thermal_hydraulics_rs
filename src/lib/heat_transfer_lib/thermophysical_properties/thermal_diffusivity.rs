use uom::si::f64::*;

use super::{Material, thermal_conductivity, density, specific_heat_capacity};


/// calculates thermal diffusivity of a material
pub fn thermal_diffusivity(material: Material, 
    temperature: ThermodynamicTemperature,
    pressure: Pressure) -> Result<DiffusionCoefficient,String> {

    let material_thermal_conductivity: ThermalConductivity = 
    thermal_conductivity::thermal_conductivity(
        material, 
        temperature, 
        pressure)?;
    
    let material_density: MassDensity = 
    density::density(
        material, 
        temperature, 
        pressure)?;
    
    let material_specific_heat_capacity: SpecificHeatCapacity = 
    specific_heat_capacity::specific_heat_capacity(
        material, 
        temperature, 
        pressure)?;

    let alpha: DiffusionCoefficient = 
    material_thermal_conductivity/ 
    material_density/ 
    material_specific_heat_capacity;

    return Ok(alpha);
}
