use uom::si::f64::*;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use uom::si::thermodynamic_temperature::kelvin;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;
use super::liquid_database::dowtherm_a::get_dowtherm_a_thermal_conductivity;

use peroxide::prelude::*;

/// returns thermal conductivity for a given material 
///
/// ```rust 
/// use uom::si::f64::*;
/// use uom::si::thermal_conductivity::watt_per_meter_kelvin;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// SolidMaterial::SteelSS304L;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// Material;
/// use thermal_hydraulics_rs::boussinesq_solver::boussinesq_thermophysical_properties::
/// thermal_conductivity::try_get_kappa_thermal_conductivity;
///
/// use uom::si::pressure::atmosphere;
///
/// let steel = Material::Solid(SteelSS304L);
/// let steel_temp = ThermodynamicTemperature::new::<kelvin>(350.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// // at 350K, we should expect thermal conductivity, 
/// // 15.58 W/(m K)
///
/// let steel_thermal_cond: ThermalConductivity = 
/// try_get_kappa_thermal_conductivity(steel, steel_temp, pressure).unwrap();
///
/// // Residuals from Graves et al. was about 3% at 350K for least 
/// // squares regression. So 2.8% error is reasonable
///
/// approx::assert_relative_eq!(
///     15.58,
///     steel_thermal_cond.value,
///     max_relative=0.028);
///
/// ``` 
#[inline]
pub fn try_get_kappa_thermal_conductivity(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<ThermalConductivity,ThermalHydraulicsLibError> {

    let thermal_conductivity: ThermalConductivity = match material {
        Material::Solid(_) => solid_thermal_conductivity(material, temperature),
        Material::Liquid(_) => liquid_thermal_conductivity(material, temperature)
    };

    return Ok(thermal_conductivity);
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
        SteelSS304L => steel_ss_304_l_ornl_thermal_conductivity(temperature),
        Copper => copper_thermal_conductivity(temperature),
    };

    return thermal_conductivity;


}

impl LiquidMaterial {
    /// returns the liquid thermal conductivity in a result enum 
    #[inline]
    pub fn try_get_thermal_conductivity(&self,
        fluid_temp: ThermodynamicTemperature,) 
        -> Result<ThermalConductivity, ThermalHydraulicsLibError>{

        let thermal_conductivity: ThermalConductivity = match self {
            DowthermA => dowtherm_a_thermal_conductivity(fluid_temp)?,
            TherminolVP1 => dowtherm_a_thermal_conductivity(fluid_temp)?
        };

        Ok(thermal_conductivity)
    }

}

impl SolidMaterial {
    /// returns the liquid thermal conductivity in a result enum 
    #[inline]
    pub fn try_get_thermal_conductivity(&self,
        solid_temp: ThermodynamicTemperature,) 
        -> Result<ThermalConductivity, ThermalHydraulicsLibError>{

        let thermal_conductivity: ThermalConductivity = match self {
        Fiberglass => fiberglass_thermal_conductivity(solid_temp) ,
        SteelSS304L => steel_ss_304_l_ornl_thermal_conductivity(solid_temp),
        Copper => copper_thermal_conductivity(solid_temp),
        };

        Ok(thermal_conductivity)
    }

}

// should the material happen to be a liquid, use this function
fn liquid_thermal_conductivity(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> ThermalConductivity {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!(
        "liquid_thermal_conductivity, use LiquidMaterial enums only")
    };

    liquid_material.try_get_thermal_conductivity(fluid_temp).unwrap()
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

    let fiberglass_thermal_conductivity_value = s.eval(
        temperature_value_kelvin);

    return ThermalConductivity::new::<watt_per_meter_kelvin>(
        fiberglass_thermal_conductivity_value);
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
#[inline]
fn steel_ss_304_l_ornl_thermal_conductivity(
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

    let copper_thermal_conductivity_value = s.
        eval(temperature_value_kelvin);

    ThermalConductivity::new::<watt_per_meter_kelvin>(
        copper_thermal_conductivity_value)

}

/// returns thermal conductivity of stainless steel 304L
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
fn _steel_304_l_spline_thermal_conductivity(
    temperature: ThermodynamicTemperature) -> ThermalConductivity {

    let temperature_value_kelvin: f64 = temperature.get::<kelvin>();
    // here we use a cubic spline to interpolate the values
    // it's a little calculation heavy, but don't really care now
    let thermal_cond_temperature_values_kelvin = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);
    let thermal_conductivity_values_watt_per_meter_kelin = c!(14.31,
        14.94, 15.58, 16.21, 16.85, 17.48, 20.02, 23.83);
    //let cp_values_watt_per_meter_kelin = c!(443.3375,
    //    457.0361, 469.4894, 480.6974, 490.66, 500.6227, 526.7746,
    //    551.6812);

    let s = CubicSpline::from_nodes(&thermal_cond_temperature_values_kelvin, 
        &thermal_conductivity_values_watt_per_meter_kelin);

    let steel_thermal_conductivity_value = s.eval(
        temperature_value_kelvin);

    return ThermalConductivity::new::<watt_per_meter_kelvin>(
        steel_thermal_conductivity_value);
}

#[inline]
fn dowtherm_a_thermal_conductivity(
    fluid_temp: ThermodynamicTemperature) -> Result<ThermalConductivity,ThermalHydraulicsLibError>{
    return get_dowtherm_a_thermal_conductivity(fluid_temp);
}

#[test]
pub fn thermal_conductivity_test_steel(){

    // we're going to test thermal conductivity for steel,
    // first at 500K for both the spline and the correlation 
    // thermal conductivity, we expect at 350K 
    // 15.58 W/(m K)

    let thermal_cond_spline = _steel_304_l_spline_thermal_conductivity(
        ThermodynamicTemperature::new::<kelvin>(350.0));

    approx::assert_relative_eq!(
        15.58,
        thermal_cond_spline.value,
        max_relative=0.001);

    // now for the Graves et al. 1991 version, from ORNL
    //

    let thermal_cond_graves_et_al_1991 = 
    steel_ss_304_l_ornl_thermal_conductivity(
        ThermodynamicTemperature::new::<kelvin>(350.0));

    // between graves and the Zou/Zweibaum version,
    // there is abut 2.8\% difference
    //
    // Residuals from Graves et al. was about 3% at 350K for least 
    // squares regression. So this is reasonable
    approx::assert_relative_eq!(
        15.58,
        thermal_cond_graves_et_al_1991.value,
        max_relative=0.028);

    // let's try now at 1000K 
    // we expect thermal thermal_conductivity to be at 23.83

    let thermal_cond_spline = 
    _steel_304_l_spline_thermal_conductivity(
        ThermodynamicTemperature::new::<kelvin>(1000.0));

    approx::assert_relative_eq!(
        23.83,
        thermal_cond_spline.value,
        max_relative=0.028);


}
