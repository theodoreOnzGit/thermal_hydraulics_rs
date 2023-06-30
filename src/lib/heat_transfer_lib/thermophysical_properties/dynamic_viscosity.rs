use uom::si::f64::*;
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::pressure::atmosphere;
use crate::fluid_mechanics_lib::therminol_component::
dowtherm_a_properties::getDowthermADensity;
use crate::fluid_mechanics_lib::therminol_component::dowtherm_a_properties::getDowthermAViscosity;
use uom::si::thermodynamic_temperature::kelvin;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;

/// returns a dynamic_viscosity given a material, temperature and pressure
///
/// example:
///
/// ```rust
/// use uom::si::f64::*;
/// use uom::si::pressure::atmosphere;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::dynamic_viscosity::dynamic_viscosity;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::LiquidMaterial::DowthermA;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::Material;
///
/// let dowtherm_a = Material::Liquid(DowthermA);
/// let temperature = ThermodynamicTemperature::new::<kelvin>(350.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// let dynamic_viscosity_result = 
/// dynamic_viscosity(dowtherm_a, temperature, pressure);
///
/// approx::assert_relative_eq!(
///     0.001237,
///     dynamic_viscosity_result.unwrap().value,
///     max_relative=0.01);
/// 
/// ```
pub fn dynamic_viscosity(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<DynamicViscosity, String> {

    match material {
        Material::Solid(_) => Err("solids don't have dynamic viscosity".to_string()),
        Material::Liquid(_) => Ok(liquid_dynamic_viscosity(material, temperature))
    }
}




// should the material happen to be a liquid, use this function
#[inline]
fn liquid_dynamic_viscosity(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> DynamicViscosity {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!("liquid_dynamic_viscosity, use LiquidMaterial enums only")
    };

    let dynamic_viscosity: DynamicViscosity = match liquid_material {
        DowthermA => dowtherm_a_dynamic_viscosity(fluid_temp),
        TherminolVP1 => dowtherm_a_dynamic_viscosity(fluid_temp)
    };

    return dynamic_viscosity;
}


#[inline]
fn dowtherm_a_dynamic_viscosity(fluid_temp: ThermodynamicTemperature) -> DynamicViscosity{
    return getDowthermAViscosity(fluid_temp);
}

#[test]
pub fn dynamic_viscosity_test_steel(){

    let dowtherm_a = Material::Liquid(DowthermA);
    let temperature = ThermodynamicTemperature::new::<kelvin>(350.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let dynamic_viscosity = dynamic_viscosity(dowtherm_a, temperature, pressure);

    approx::assert_relative_eq!(
        0.001237,
        dynamic_viscosity.unwrap().value,
        max_relative=0.01);
}
