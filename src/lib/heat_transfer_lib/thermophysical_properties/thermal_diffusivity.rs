use uom::si::f64::*;

use super::{Material, thermal_conductivity, density, specific_heat_capacity};


/// calculates thermal diffusivity of a material
/// ```rust
/// use uom::si::f64::*;
/// use uom::si::pressure::atmosphere;
/// use uom::si::thermodynamic_temperature::degree_celsius;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::thermal_diffusivity::try_get_alpha_thermal_diffusivity;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::thermal_conductivity::try_get_kappa_thermal_conductivity;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::SolidMaterial::SteelSS304L;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::Material;
///
/// let steel = Material::Solid(SteelSS304L);
/// let temperature = ThermodynamicTemperature::new
/// ::<degree_celsius>(80.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// let thermal_diffusivity_result = try_get_alpha_thermal_diffusivity(
/// steel, temperature, pressure).unwrap();
///  
/// // thermal diffusivity of ss304L is approx 4.13e-6 m^2/s
/// approx::assert_relative_eq!(
/// thermal_diffusivity_result.value,
/// 4.13e-6,
/// epsilon = 0.001);
///
/// // conductivity is approx 15.62 W/(m K)
/// let steel_thermal_cond: ThermalConductivity = 
/// try_get_kappa_thermal_conductivity(steel, temperature, pressure).unwrap();
/// 
/// approx::assert_relative_eq!(
/// steel_thermal_cond.value,
/// 15.62,
/// epsilon = 0.001);
/// ```
#[inline]
pub fn try_get_alpha_thermal_diffusivity(material: Material, 
    temperature: ThermodynamicTemperature,
    pressure: Pressure) -> Result<DiffusionCoefficient,String> {

    let material_thermal_conductivity: ThermalConductivity = 
    thermal_conductivity::try_get_kappa_thermal_conductivity(
        material, 
        temperature, 
        pressure)?;
    
    let material_density: MassDensity = 
    density::try_get_rho(
        material, 
        temperature, 
        pressure)?;
    
    let material_specific_heat_capacity: SpecificHeatCapacity = 
    specific_heat_capacity::try_get_cp(
        material, 
        temperature, 
        pressure)?;

    let alpha: DiffusionCoefficient = 
    material_thermal_conductivity/ 
    material_density/ 
    material_specific_heat_capacity;

    return Ok(alpha);
}

