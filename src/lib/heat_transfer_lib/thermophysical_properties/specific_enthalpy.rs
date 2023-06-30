use uom::si::available_energy::joule_per_gram;
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
/// use uom::si::specific_heat_capacity::{joule_per_kilogram_kelvin,
/// joule_per_gram_degree_celsius};
/// use uom::si::thermodynamic_temperature::kelvin;
/// use uom::si::temperature_interval::degree_celsius;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// SolidMaterial::{SteelSS304L,Copper};
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// Material;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// specific_enthalpy::specific_enthalpy;
///
/// use uom::si::pressure::atmosphere;
///
/// let steel = Material::Solid(SteelSS304L);
/// let steel_temp = ThermodynamicTemperature::new::<kelvin>(273.15);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// // enthalpy should be zero at 273.15 K
///
/// let steel_enthalpy_273_15_kelvin = 
/// specific_enthalpy(steel, steel_temp, pressure);
///
/// approx::assert_relative_eq!(
///     0.0,
///     steel_enthalpy_273_15_kelvin.unwrap().value,
///     max_relative=0.045);
/// 
/// // we can also calculate enthalpy change of copper 
/// // from 375K to 425K
/// let test_temperature_1 = ThermodynamicTemperature::new::
/// <kelvin>(375.0);
/// let test_temperature_2 = ThermodynamicTemperature::new::
/// <kelvin>(425.0);
///
/// let copper = Material::Solid(Copper);
///
/// let copper_enthalpy_change = 
/// specific_enthalpy(copper, test_temperature_2, pressure).unwrap()
/// - specific_enthalpy(copper, test_temperature_1, pressure).unwrap();
///
/// // http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/sphtt.html
/// // https://www.engineeringtoolbox.com/specific-heat-metals-d_152.html
/// // copper at 20C has heat capacity of 
/// // 0.386 J/(g K)
/// // going to use this to estimate a ballpark figure to find enthalpy 
/// // h = cp(T2 - T1)
/// 
/// // we can't usually subtract thermodynamic temperatures from each 
/// // other, we need a termpature interval
/// // 
///
/// let cp_copper_20_c = 
/// SpecificHeatCapacity::new::<joule_per_gram_degree_celsius>(0.386);
/// 
/// let temperature_difference = 
/// TemperatureInterval::new::<degree_celsius>(
/// test_temperature_2.value - test_temperature_1.value);
///
/// let specific_enthalpy_ballpark = 
/// cp_copper_20_c * temperature_difference;
/// 
/// // the ballpark value is 19300 J/kg
/// approx::assert_relative_eq!(
///     specific_enthalpy_ballpark.value,
///     19300.0,
///     max_relative=0.0001);
///
/// // it's less than 4% different from the ballpark value
/// // This means the copper enthalpy change should be quite reasonable
///
/// approx::assert_relative_eq!(
///     specific_enthalpy_ballpark.value,
///     copper_enthalpy_change.value,
///     max_relative=0.04);
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
//
// probably should have a temperature range checker in 
// future
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

///
/// Graves, R. S., Kollie, T. G., 
/// McElroy, D. L., & Gilchrist, K. E. (1991). The 
/// thermal conductivity of AISI 304L stainless steel. 
/// International journal of thermophysics, 12, 409-415. 
///
/// data taken from ORNL
///
/// It's only good for range of 300K to 700K
///
/// However, I analytically integrated it with wolfram alpha
#[inline]
fn steel_ss_304_l_ornl_specific_enthalpy(
    temperature: ThermodynamicTemperature) -> AvailableEnergy {

    // first I define a function for specific enthalpy between two 
    // temperatures in kelvin
    fn definite_integral_specific_enthalpy(
        temp_1: ThermodynamicTemperature,
        temp_2: ThermodynamicTemperature) -> AvailableEnergy {

        // the integration constant is assumed to be zero 

        let temp_1_value_kelvin = temp_1.get::<kelvin>();
        let temp_2_value_kelvin = temp_2.get::<kelvin>();

        let enthalpy_value_joule_per_gram_per_kelvin_temp_1 = 
        1.73333e-8 * f64::powf(temp_1_value_kelvin,3.0) 
        + 0.000085 * f64::powf(temp_1_value_kelvin, 2.0)
        + 0.4267 * temp_1_value_kelvin;

        let enthalpy_value_joule_per_gram_per_kelvin_temp_2 = 
        1.73333e-8 * f64::powf(temp_2_value_kelvin,3.0) 
        + 0.000085 * f64::powf(temp_2_value_kelvin, 2.0)
        + 0.4267 * temp_2_value_kelvin;

        let enthalpy_difference_joule_per_gram = 
        enthalpy_value_joule_per_gram_per_kelvin_temp_2 
        - enthalpy_value_joule_per_gram_per_kelvin_temp_1;

        AvailableEnergy::new::<joule_per_gram>(
            enthalpy_difference_joule_per_gram)
    }

    // reference temperature is zero degrees c, 
    // enthalpy is zero j/kg at that point
    let refernce_temperature = ThermodynamicTemperature::new::
    <kelvin>(273.15);

    let steel_enthalpy = definite_integral_specific_enthalpy(
        refernce_temperature, temperature);

    steel_enthalpy
}

#[test]
pub fn specific_enthalpy_test_steel_ornl(){
    // let's test specifc enthalpy at 350K 

    let test_temperature = ThermodynamicTemperature::new::
    <kelvin>(350.0);

    // wolfram gives an enthalpy (assuming enthalpy is zero at zero 
    // degrees C, 273.15 K)
    // this is done using the Graves et al. 1991 version for cp
    //  37.2524 j/g
    let wolfram_enthalpy_value_joule_per_kg = 37.2524*1000.0;

    let enthalpy_analytical_ornl = 
    steel_ss_304_l_ornl_specific_enthalpy(test_temperature);

    approx::assert_relative_eq!(
        wolfram_enthalpy_value_joule_per_kg,
        enthalpy_analytical_ornl.value,
        max_relative=0.0001);

    
}

/// here is a test for comparing ornl and nico zweibaum's value 
/// at 375 to 425 kelvin
///
/// the cp correlation was from 300 to 700 Kelvin, so using 273.15 
/// as zero enthalpy is technically outside the range.
///
/// Despite this, it should still be able to calculate enthalpy 
/// change from 375 K to 425K
#[test]
pub fn specific_enthalpy_test_steel_ornl_and_zweibaum_spline(){
    // let's test specifc enthalpy at 350K 

    let test_temperature_1 = ThermodynamicTemperature::new::
    <kelvin>(375.0);
    let test_temperature_2 = ThermodynamicTemperature::new::
    <kelvin>(425.0);

    // wolfram gives an enthalpy (assuming enthalpy is zero at zero 
    // degrees C, 273.15 K)
    // this is done using the Graves et al. 1991 version for cp
    // 25.1515 j/g for 375 to 425 K
    let wolfram_enthalpy_value_joule_per_kg = 25.1515*1000.0;

    let enthalpy_analytical_ornl = 
    steel_ss_304_l_ornl_specific_enthalpy(test_temperature_2)
    - steel_ss_304_l_ornl_specific_enthalpy(test_temperature_1);

    approx::assert_relative_eq!(
        wolfram_enthalpy_value_joule_per_kg,
        enthalpy_analytical_ornl.value,
        max_relative=0.0001);

    // now let's test the spline version 
    //
    let enthalpy_spline_zweibaum = 
    steel_304_l_spline_specific_enthalpy(test_temperature_2)
    - steel_304_l_spline_specific_enthalpy(test_temperature_1);

    // there is about a 4.5% difference between the ornl value 
    // and the spline value
    // doesn't seem too bad honestly.
    //
    // Of course, one can do uncertainty propagation in order to 
    // find out the degree of change, but I won't do that for now.
    //
    // otherwise, spline should work quite okay
    approx::assert_relative_eq!(
        enthalpy_analytical_ornl.value,
        enthalpy_spline_zweibaum.value,
        max_relative=0.045);

    


    
}

