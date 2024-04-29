use uom::si::f64::*;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use uom::si::thermodynamic_temperature::kelvin;

use super::range_check;
use super::solid_database::ss_304_l::steel_304_l_libreoffice_spline_specific_heat_capacity_ciet_zweibaum;
use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;
use super::liquid_database::dowtherm_a::get_dowtherm_a_constant_pressure_specific_heat_capacity;

use peroxide::prelude::*;

/// returns cp for a given material 
///
/// ```rust 
/// use uom::si::f64::*;
/// use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// SolidMaterial::SteelSS304L;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// Material;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// specific_heat_capacity::try_get_cp;
///
/// use uom::si::pressure::atmosphere;
///
/// let steel = Material::Solid(SteelSS304L);
/// let steel_temp = ThermodynamicTemperature::new::<kelvin>(350.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// // at 350K, we should expect thermal conductivity, 
/// // 470 W/(m K)
///
/// let steel_thermal_cond: SpecificHeatCapacity = 
/// try_get_cp(steel, steel_temp, pressure).unwrap();
///
///
/// approx::assert_relative_eq!(
///     470.0,
///     steel_thermal_cond.value,
///     max_relative=0.035);
///
/// ``` 
#[inline]
pub fn try_get_cp(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<SpecificHeatCapacity, ThermalHydraulicsLibError> {

    let specific_heat_capacity: SpecificHeatCapacity = match material {
        Material::Solid(_) => solid_specific_heat_capacity(material, temperature)?,
        Material::Liquid(_) => liquid_specific_heat_capacity(material, temperature)?
    };

    return Ok(specific_heat_capacity);
}

// should the material happen to be a solid, use this function
fn solid_specific_heat_capacity(material: Material,
    temperature: ThermodynamicTemperature) -> Result<SpecificHeatCapacity, ThermalHydraulicsLibError>{
    
    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Liquid(_) => panic!("solid_specific_heat_capacity, use SolidMaterial enums only")
    };

    let specific_heat_capacity: SpecificHeatCapacity = match solid_material {
        Fiberglass => fiberglass_specific_heat_capacity(temperature) ,
        SteelSS304L => steel_304_l_libreoffice_spline_specific_heat_capacity_ciet_zweibaum(temperature)?,
        Copper => copper_specific_heat_capacity(temperature)?,
    };

    return Ok(specific_heat_capacity);


}

// should the material happen to be a liquid, use this function
fn liquid_specific_heat_capacity(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> Result<SpecificHeatCapacity,
ThermalHydraulicsLibError>{

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!(
        "liquid_specific_heat_capacity, use LiquidMaterial enums only")
    };

    let specific_heat_capacity: SpecificHeatCapacity = match liquid_material {
        DowthermA => dowtherm_a_specific_heat_capacity(fluid_temp)?,
        TherminolVP1 => dowtherm_a_specific_heat_capacity(fluid_temp)?
    };

    return Ok(specific_heat_capacity);
}

/// returns thermal conductivity of fiberglass
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn fiberglass_specific_heat_capacity(
    _temperature: ThermodynamicTemperature) -> SpecificHeatCapacity {

    return SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
        844.0);
}


/// returns thermal conductivity of copper
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn copper_specific_heat_capacity(
    temperature: ThermodynamicTemperature) -> Result<SpecificHeatCapacity,ThermalHydraulicsLibError> {

    range_check(
        &Material::Solid(SolidMaterial::Copper),
        temperature, 
        ThermodynamicTemperature::new::<kelvin>(1000.0), 
        ThermodynamicTemperature::new::<kelvin>(200.0))?;

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy
    //
    // and actually, rebuilding the spline is quite problematic
    // we need to build it ONCE and read from it
    //
    let thermal_cond_temperature_values_kelvin = c!(200.0, 
        250.0, 300.0, 350.0, 
        400.0, 500.0, 1000.0);
    let specific_heat_capacity_values_joule_per_kilogram_kelvin = c!(355.7047,
        373.6018, 384.7875, 392.6174,
        398.2103, 407.1588, 417.2260);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &specific_heat_capacity_values_joule_per_kilogram_kelvin);

    let copper_specific_heat_capacity_value = s.
        eval(temperature_value_kelvin);

    Ok(SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
        copper_specific_heat_capacity_value))

}


#[inline]
fn dowtherm_a_specific_heat_capacity(
    fluid_temp: ThermodynamicTemperature) -> Result<SpecificHeatCapacity
,ThermalHydraulicsLibError>{
    get_dowtherm_a_constant_pressure_specific_heat_capacity(fluid_temp)
}

