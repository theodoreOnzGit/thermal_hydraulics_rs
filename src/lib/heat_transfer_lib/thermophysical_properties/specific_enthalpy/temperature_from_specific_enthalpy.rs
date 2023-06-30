use uom::si::available_energy::joule_per_gram;
use uom::si::f64::*;
use uom::si::available_energy::joule_per_kilogram;
use crate::fluid_mechanics_lib::therminol_component::
dowtherm_a_properties::getDowthermAEnthalpy;
use uom::si::thermodynamic_temperature::kelvin;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;
use uom::si::pressure::atmosphere;
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy
::specific_enthalpy;

use peroxide::prelude::*;

// should the material happen to be a solid, use this function
// to extract temperature from enthalpy
//
// probably should have a temperature range checker in 
// future
//
// pub(in crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy) 
// here only makes it accessible to the 
// specific_enthalpy/mod.rs 
// nothing else
pub(in crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy) 
fn get_solid_temperature_from_specific_enthalpy(material: Material,
    h_material: AvailableEnergy) -> ThermodynamicTemperature {
    
    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Liquid(_) => panic!("solid_specific_enthalpy, use SolidMaterial enums only")
    };

    let material_temperature: ThermodynamicTemperature = 
    match solid_material {
        Fiberglass => todo!(),
        SteelSS304L => 
            steel_304_l_spline_temperature_from_specific_enthalpy(
                h_material),
        Copper => todo!(),
    };

    return material_temperature;


}

/// returns temperature of stainless steel 304L 
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// The algorithm is to make a spline of enthalpy and temperature
///
/// Note: h in this function represents specific enthalpy
#[inline]
fn steel_304_l_spline_temperature_from_specific_enthalpy(
    h_steel: AvailableEnergy) -> ThermodynamicTemperature {

    let temperature_values_kelvin: Vec<f64>
    = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);

    let specific_heat_capacity_values_joule_per_kilogram_kelvin = c!(443.3375,
        457.0361, 469.4894, 480.6974, 490.66, 500.6227, 526.7746,
        551.6812);

    let h_spline_integral = 
    CubicSpline::from_nodes(&temperature_values_kelvin, 
        &specific_heat_capacity_values_joule_per_kilogram_kelvin)
        .integral();

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let zero_celsius_vec = vec![273.15; temperature_vec_len];

    // now the integrals are indefinite, hence I need to ensure 
    // the enthalpies are taken with reference to the common 
    // reference point of 0 degrees C or 273.15 K
    let integral_spline_enthalpy_vector = 
    h_spline_integral.eval_vec(&temperature_values_kelvin);
    let reference_enthalpy_vector =
    h_spline_integral.eval_vec(&zero_celsius_vec);

    // this is a concise but rather complicated way of 
    // subtraction
    let enthalpy_vector: Vec<f64> = integral_spline_enthalpy_vector
        .iter()
        .zip(reference_enthalpy_vector)
        .map(|(elem_a, elem_b)| elem_a - elem_b)
        .collect();
    
    // now I have my enthalpy vector, i can do an inverted spline 
    // to have enthalpy given in as an input, and temperature received
    // as an output

    let enthalpy_to_temperature_spline = 
    CubicSpline::from_nodes(&enthalpy_vector,
    &temperature_values_kelvin);

    // now let's get our enthalpy in joules_per_kg
    let h_steel_joules_per_kg = h_steel.get::<joule_per_kilogram>();

    let temperature_from_enthalpy_kelvin = 
    enthalpy_to_temperature_spline.eval(h_steel_joules_per_kg);

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}

#[test]
pub fn temperature_from_enthalpy_test_spline(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let steel = Material::Solid(SteelSS304L);
    let steel_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        steel,steel_temp,pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    steel_304_l_spline_temperature_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically off by about 21K,
    // the temperature was i think 354K
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=0.0);

    // looks like the spline is good as an initial guess, but not 
    // really for finding the answer

}
