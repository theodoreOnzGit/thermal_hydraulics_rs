use uom::si::f64::*;
use uom::si::available_energy::joule_per_kilogram;
use crate::fluid_mechanics_lib::therminol_component::dowtherm_a_properties::getDowthermAEnthalpy;
use uom::si::thermodynamic_temperature::kelvin;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;

use peroxide::prelude::*;

/// returns specific enthaply for a given material 
/// specific_enthalpy is defined as 0 for 0 degree_celsius
/// for any material, that is 273.15 K
///
/// ```rust 
/// use uom::si::f64::*;
/// use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// SolidMaterial::SteelSS304L;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// Material;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// specific_heat_capacity::specific_heat_capacity;
///
/// use uom::si::pressure::atmosphere;
///
/// let steel = Material::Solid(SteelSS304L);
/// let steel_temp = ThermodynamicTemperature::new::<kelvin>(350.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
/// todo!("write doctest for specific enthalpy");
///
/// ``` 
pub fn specific_enthalpy(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<AvailableEnergy, String> {

    let specific_enthalpy: AvailableEnergy = match material {
        Material::Solid(_) => solid_specific_enthalpy(material, temperature),
        Material::Liquid(_) => liquid_specific_enthalpy(material, temperature)
    };

    return Ok(specific_enthalpy);
}

// should the material happen to be a solid, use this function
fn solid_specific_enthalpy(material: Material,
    temperature: ThermodynamicTemperature) -> AvailableEnergy {
    
    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Liquid(_) => panic!("solid_specific_enthalpy, use SolidMaterial enums only")
    };

    let specific_enthalpy: AvailableEnergy = match solid_material {
        Fiberglass => fiberglass_specific_enthalpy(temperature) ,
        SteelSS304L => steel_304_l_spline_specific_enthalpy(temperature),
        Copper => copper_specific_enthalpy(temperature),
    };

    return specific_enthalpy;


}

// should the material happen to be a liquid, use this function
fn liquid_specific_enthalpy(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> AvailableEnergy {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!(
        "liquid_specific_enthalpy, use LiquidMaterial enums only")
    };

    let specific_enthalpy: AvailableEnergy = match liquid_material {
        DowthermA => dowtherm_a_specific_enthalpy(fluid_temp),
        TherminolVP1 => dowtherm_a_specific_enthalpy(fluid_temp)
    };

    return specific_enthalpy;
}

/// returns thermal conductivity of fiberglass
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// specific enthalpy at 273.15 K = 0
#[inline]
fn fiberglass_specific_enthalpy(
    temperature: ThermodynamicTemperature) -> AvailableEnergy {

    let temp_value_kelvin = temperature.get::<kelvin>();
    // fiberglass cp in J/kg.K
    fn fiberglass_cp_float(temp_value_kelvin: f64) -> f64 {

        let temperature = ThermodynamicTemperature::
            new::<kelvin>(temp_value_kelvin);
        let fiberglass_specific_heat_capacity_value_j_per_kg_k: f64 = 
        fiberglass_specific_enthalpy(temperature).value;

        fiberglass_specific_heat_capacity_value_j_per_kg_k
    }

    let specific_enthalpy_value_j_per_kg = 
    integrate(fiberglass_cp_float, (273.15,temp_value_kelvin));

    return AvailableEnergy::new::<joule_per_kilogram>(
        specific_enthalpy_value_j_per_kg);
}


/// returns thermal conductivity of copper
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn copper_specific_enthalpy(
    temperature: ThermodynamicTemperature) -> AvailableEnergy {

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy, but don't really care now
    let thermal_cond_temperature_values_kelvin = c!(200.0, 
        250.0, 300.0, 350.0, 
        400.0, 500.0, 1000.0);
    let specific_heat_capacity_values_joule_per_kilogram_kelvin = c!(355.7047,
        373.6018, 384.7875, 392.6174,
        398.2103, 407.1588, 417.2260);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &specific_heat_capacity_values_joule_per_kilogram_kelvin);

    let copper_specific_enthalpy_value = s.
        integrate((273.15,temperature_value_kelvin));

    return AvailableEnergy::new::<joule_per_kilogram>(
        copper_specific_enthalpy_value);

}

/// returns thermal conductivity of stainless steel 304L
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn steel_304_l_spline_specific_enthalpy(
    temperature: ThermodynamicTemperature) -> AvailableEnergy {


    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy, but don't really care now
    let thermal_cond_temperature_values_kelvin = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);
    let specific_heat_capacity_values_joule_per_kilogram_kelvin = c!(443.3375,
        457.0361, 469.4894, 480.6974, 490.66, 500.6227, 526.7746,
        551.6812);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &specific_heat_capacity_values_joule_per_kilogram_kelvin);

    let steel_specific_enthalpy_value = s.integrate(
        (273.15,temperature_value_kelvin));

    return AvailableEnergy::new::<joule_per_kilogram>(
        steel_specific_enthalpy_value);
}

#[inline]
fn dowtherm_a_specific_enthalpy(
    fluid_temp: ThermodynamicTemperature) -> AvailableEnergy{
    return getDowthermAEnthalpy(fluid_temp);
}

#[test]
pub fn specific_enthalpy_test_steel(){
    todo!("need to write test for specific enthalpy")
}
