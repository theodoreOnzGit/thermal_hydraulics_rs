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

// for spline method
use peroxide::prelude::*;

// for root finding with brent
extern crate roots;
use roots::Roots;
use roots::find_root_brent;
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
            steel_304_l_spline_temp_attempt_3_from_specific_enthalpy(
                h_material),
        Copper => 
            copper_spline_temp_attempt_2_from_specific_enthalpy(
                h_material),
    };

    return material_temperature;


}
#[inline]
fn copper_spline_temp_attempt_2_from_specific_enthalpy(
    h_copper: AvailableEnergy) -> ThermodynamicTemperature {

    // the idea is basically to evaluate enthalpy at the 
    // following temperatures
    let temperature_values_kelvin: Vec<f64>
    = c!(200.0 ,250.0, 300.0, 350.0, 
        400.0, 500.0, 1000.0);

    // and then use that to formulate a spline,
    // with the spline, i'll evaluate enthalpy from temperature
    // within pretty much one iteration. However, it is spline 
    // construction which may take a little long. 
    //
    // However, the number of iterations per calculation is fixed
    //
    // I won't optimise it now just yet

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let mut enthalpy_vector = vec![0.0; temperature_vec_len];

    for index_i in 0..temperature_vec_len {

        // first, evaluate the enthalpy at temperature values 
        let temperature_value = temperature_values_kelvin[index_i];

        //next let's evaluate the specific enthalpy of copper 
        let copper = Material::Solid(Copper);
        let copper_temp = ThermodynamicTemperature::new::<kelvin>(
            temperature_value);
        let pressure = Pressure::new::<atmosphere>(1.0);

        let copper_enthalpy_result = specific_enthalpy(copper, 
            copper_temp, pressure);

        let copper_enthalpy_value = match copper_enthalpy_result {
            Ok(copper_enthalpy) => copper_enthalpy.value,
            // i can of course unwrap the result,
            // but i want to leave it more explicit in case 
            // i wish to manually handle the error
            Err(error_msg) => panic!("{}",error_msg),
        };

        // once i evalute the enthalpy value, pass it on to the vector

        enthalpy_vector[index_i] = copper_enthalpy_value;

    }


    // now I have my enthalpy vector, i can do an inverted spline 
    // to have enthalpy given in as an input, and temperature received
    // as an output

    let enthalpy_to_temperature_spline = 
    CubicSpline::from_nodes(&enthalpy_vector,
    &temperature_values_kelvin);

    // now let's get our enthalpy in joules_per_kg
    let h_copper_joules_per_kg = h_copper.get::<joule_per_kilogram>();

    let temperature_from_enthalpy_kelvin = 
    enthalpy_to_temperature_spline.eval(h_copper_joules_per_kg);

    // now, the copper enthalpy will not be quite near 
    // enough, but it is very close. We can bracket 
    // the root 


    let enthalpy_root = |temp_degrees_c_value : f64| -> f64 {
        let lhs_value = h_copper.get::<joule_per_kilogram>();


        let copper = Material::Solid(Copper);
        let copper_temp = ThermodynamicTemperature::new::
            <kelvin>(temp_degrees_c_value) ;
        let pressure = Pressure::new::<atmosphere>(1.0);

        let rhs = specific_enthalpy(copper, 
            copper_temp, pressure);

        let rhs_value = match rhs {
            Ok(enthalpy_val) => enthalpy_val.get::<joule_per_kilogram>(),
                // fall back to guess value, 
            Err(error_msg) => panic!("{}",error_msg),
        };

        return lhs_value-rhs_value;
    };

    let brent_error_bound: f64 = 30.0;

    let upper_limit: f64 = temperature_from_enthalpy_kelvin +
        brent_error_bound;

    let lower_limit : f64 = temperature_from_enthalpy_kelvin -
        brent_error_bound;

    let fluid_temperature_degrees_c_result
    = find_root_brent(upper_limit,
        lower_limit,
        enthalpy_root,
        &mut 1e-8_f64
    );

    let temperature_from_enthalpy_kelvin = 
    fluid_temperature_degrees_c_result.unwrap();

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}


#[test]
pub fn copper_temperature_from_enthalpy_test_spline_2(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let copper = Material::Solid(Copper);
    let copper_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        copper,copper_temp,pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    copper_spline_temp_attempt_2_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically by about 5K, which is 
    // not within measurement error, probably have to do more work
    // what this means is that accuracy is sacrificed
    // for speed, sometimes too much accuracy
    //
    // for enthalpy, we probably want to have it as accurate 
    // as possible so that energy doesn't appear from nowhere 
    // and disappear from calculation
    //
    // I would note though, that the spline method does 
    // give a pretty good initial guess of where the temperature 
    // ought to be, so perhaps the iterative method can be used 
    // for the last few iterations to convergence
    // we could use brent dekker method
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=0.005);


}


/// returns temperature of stainless steel 304L 
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// The algorithm is to make a spline of enthalpy and temperature
/// and then use the enthalpy to obtain a temperature
///
/// Note: h in this function represents specific enthalpy
///
/// This uses the spline methodology as an initial guess
/// but then uses brent dekker to finish it off
#[inline]
fn steel_304_l_spline_temp_attempt_3_from_specific_enthalpy(
    h_steel: AvailableEnergy) -> ThermodynamicTemperature {

    // the idea is basically to evaluate enthalpy at the 
    // following temperatures
    let temperature_values_kelvin: Vec<f64>
    = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);

    // and then use that to formulate a spline,
    // with the spline, i'll evaluate enthalpy from temperature
    // within pretty much one iteration. However, it is spline 
    // construction which may take a little long. 
    //
    // However, the number of iterations per calculation is fixed
    //
    // I won't optimise it now just yet

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let mut enthalpy_vector = vec![0.0; temperature_vec_len];

    for index_i in 0..temperature_vec_len {

        // first, evaluate the enthalpy at temperature values 
        let temperature_value = temperature_values_kelvin[index_i];

        //next let's evaluate the specific enthalpy of steel 
        let steel = Material::Solid(SteelSS304L);
        let steel_temp = ThermodynamicTemperature::new::<kelvin>(
            temperature_value);
        let pressure = Pressure::new::<atmosphere>(1.0);

        let steel_enthalpy_result = specific_enthalpy(steel, 
            steel_temp, pressure);

        let steel_enthalpy_value = match steel_enthalpy_result {
            Ok(steel_enthalpy) => steel_enthalpy.value,
            // i can of course unwrap the result,
            // but i want to leave it more explicit in case 
            // i wish to manually handle the error
            Err(error_msg) => panic!("{}",error_msg),
        };

        // once i evalute the enthalpy value, pass it on to the vector

        enthalpy_vector[index_i] = steel_enthalpy_value;

    }


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

    let enthalpy_root = |temp_degrees_c_value : f64| -> f64 {
        let lhs_value = h_steel.get::<joule_per_kilogram>();


        let copper = Material::Solid(SteelSS304L);
        let copper_temp = ThermodynamicTemperature::new::
            <kelvin>(temp_degrees_c_value) ;
        let pressure = Pressure::new::<atmosphere>(1.0);

        let rhs = specific_enthalpy(copper, 
            copper_temp, pressure);

        let rhs_value = match rhs {
            Ok(enthalpy_val) => enthalpy_val.get::<joule_per_kilogram>(),
                // fall back to guess value, 
            Err(error_msg) => panic!("{}",error_msg),
        };

        return lhs_value-rhs_value;
    };

    // brent error bounds can be smaller for this steel spline, 
    // it's a lot more accurate than for copper
    let brent_error_bound: f64 = 1.0;

    let upper_limit: f64 = temperature_from_enthalpy_kelvin +
        brent_error_bound;

    let lower_limit : f64 = temperature_from_enthalpy_kelvin -
        brent_error_bound;

    let fluid_temperature_degrees_c_result
    = find_root_brent(upper_limit,
        lower_limit,
        enthalpy_root,
        &mut 1e-8_f64
    );

    let temperature_from_enthalpy_kelvin = 
    fluid_temperature_degrees_c_result.unwrap();

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}
#[test]
pub fn steel_temperature_from_enthalpy_test_spline_3(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let steel = Material::Solid(SteelSS304L);
    let steel_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        steel,steel_temp,atmospheric_pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    steel_304_l_spline_temp_attempt_3_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically off by less than 0.05K, which is 
    // within measurement error!
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=0.00005);

    // let's test from 325, 425, 525, 625, 725, 825, 925
    let temperature_vec_kelvin: Vec<f64> = 
    vec![325.0, 425.0, 525.0, 625.0, 725.0, 825.0, 925.0];

    for temperature_val_kelvin in temperature_vec_kelvin.iter() {

        let steel_temp = ThermodynamicTemperature::new::<kelvin>(
            *temperature_val_kelvin);

        let enthalpy_spline_zweibaum = specific_enthalpy(steel, 
            steel_temp, atmospheric_pressure).unwrap();

        let temperature_from_enthalpy_test = 
        steel_304_l_spline_temp_attempt_3_from_specific_enthalpy(
            enthalpy_spline_zweibaum);

        // we are basically off by less than 0.5K, which is 
        // within measurement error!
        approx::assert_abs_diff_eq!(
            temperature_from_enthalpy_test.get::<kelvin>(),
            *temperature_val_kelvin,
            epsilon=0.00005);
    };


}

