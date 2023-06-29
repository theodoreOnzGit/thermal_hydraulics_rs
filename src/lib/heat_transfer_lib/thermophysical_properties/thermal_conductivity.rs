use uom::si::f64::*;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use uom::si::pressure::atmosphere;
use crate::fluid_mechanics_lib::therminol_component::
dowtherm_a_properties::getDowthermADensity;
use crate::fluid_mechanics_lib::therminol_component::dowtherm_a_properties::getDowthermAThermalConductivity;
use uom::si::thermodynamic_temperature::kelvin;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;

use peroxide::prelude::*;

pub fn thermal_conductivity(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> ThermalConductivity {

    let thermal_conductivity: ThermalConductivity = match material {
        Material::Solid(_) => solid_thermal_conductivity(material, temperature),
        Material::Liquid(_) => liquid_thermal_conductivity(material, temperature)
    };

    return thermal_conductivity;
}

// should the material happen to be a solid, use this function
fn solid_thermal_conductivity(material: Material,
    temperature: ThermodynamicTemperature) -> ThermalConductivity{
    
    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Liquid(_) => panic!("solid_thermal_conductivity, use SolidMaterial enums only")
    };

    let thermal_conductivity: ThermalConductivity = match solid_material {
        Fiberglass => fiberglass_thermal_conductivity(temperature) ,
        SteelSS304L => steel_ss_304_l_thermal_conductivity(temperature),
        Copper => copper_thermal_conductivity(temperature),
    };

    return thermal_conductivity;


}

// should the material happen to be a liquid, use this function
fn liquid_thermal_conductivity(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> ThermalConductivity {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!("liquid_thermal_conductivity, use LiquidMaterial enums only")
    };

    let thermal_conductivity: ThermalConductivity = match liquid_material {
        DowthermA => dowtherm_a_thermal_conductivity(fluid_temp),
        TherminolVP1 => dowtherm_a_thermal_conductivity(fluid_temp)
    };

    return thermal_conductivity;
}

/// returns thermal conductivity of fiberglass
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn fiberglass_thermal_conductivity(
    temperature: ThermodynamicTemperature) -> ThermalConductivity {

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy, but don't really care now
    let thermal_cond_temperature_values_kelvin = c!(250.0, 293.15, 350.0, 
        400.0, 500.0, 600.0);
    let thermal_conductivity_values_watt_per_meter_kelin = c!(0.028616,
        0.033060, 0.038916, 0.044066, 0.054366, 0.064666);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &thermal_conductivity_values_watt_per_meter_kelin);

    let copper_thermal_conductivity_value = s.polynomial_at(
        temperature_value_kelvin);

    return ThermalConductivity::new::<watt_per_meter_kelvin>(
        temperature_value_kelvin);
}

///
/// Graves, R. S., Kollie, T. G., 
/// McElroy, D. L., & Gilchrist, K. E. (1991). The 
/// thermal conductivity of AISI 304L stainless steel. 
/// International journal of thermophysics, 12, 409-415. 
#[inline]
fn steel_ss_304_l_thermal_conductivity(
    temperature: ThermodynamicTemperature) -> ThermalConductivity {

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    let thermal_conductivity_val = 7.9318 
    + 0.023051 * temperature_value_kelvin
    - 6.4166 * f64::powf(10.0, -6.0);

    ThermalConductivity::new::<watt_per_meter_kelvin>(
        thermal_conductivity_val)
}

/// returns thermal conductivity of copper
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn copper_thermal_conductivity(
    temperature: ThermodynamicTemperature) -> ThermalConductivity {

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy, but don't really care now
    let thermal_cond_temperature_values_kelvin = c!(250.0, 300.0, 350.0, 
        400.0, 500.0, 1000.0);
    let thermal_conductivity_values_watt_per_meter_kelin = c!(406.0,
        401.0, 369.0, 393.0, 386.0, 352.0);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &thermal_conductivity_values_watt_per_meter_kelin);

    let copper_thermal_conductivity_value = s.polynomial_at(
        temperature_value_kelvin);

    return ThermalConductivity::new::<watt_per_meter_kelvin>(
        temperature_value_kelvin);
}

#[inline]
fn dowtherm_a_thermal_conductivity(
    fluid_temp: ThermodynamicTemperature) -> ThermalConductivity{
    return getDowthermAThermalConductivity(fluid_temp);
}

#[test]
pub fn thermal_conductivity_test_steel(){

    todo!();
}
