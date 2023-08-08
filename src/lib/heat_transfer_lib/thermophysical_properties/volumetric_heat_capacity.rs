use uom::si::f64::*;

use super::{Material, density, specific_heat_capacity};


/// calculates volumetric_heat_capacity of a material
pub fn rho_cp(material: Material, 
    temperature: ThermodynamicTemperature,
    pressure: Pressure) -> Result<VolumetricHeatCapacity,String> {

    
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

    let rho_cp: VolumetricHeatCapacity = 
    material_density*
    material_specific_heat_capacity;

    return Ok(rho_cp);
}

