use peroxide::fuga::{CubicSpline, Spline};
use uom::si::f64::*;
use uom::si::length::micrometer;
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use uom::si::thermodynamic_temperature::kelvin;

/// density ranges not quite given in original text 
/// Zou, Ling, Rui Hu, and Anne Charpentier. SAM code validation 
/// using the compact integral effects test (CIET) experimental data. 
/// No. ANL/NSE-19/11. 
/// Argonne National Lab.(ANL), Argonne, IL (United States), 2019.
#[inline]
pub fn copper_density() -> Result<MassDensity,ThermalHydraulicsLibError> {
    return Ok(MassDensity::new::<kilogram_per_cubic_meter>(8940.0));
}

/// Arenales, M. R. M., Kumar, S., 
/// Kuo, L. S., & Chen, P. H. (2020). 
/// Surface roughness variation effects on copper tubes in 
/// pool boiling of water. International Journal of 
/// Heat and Mass Transfer, 151, 119399.
pub fn copper_surf_roughness() -> Length {
    Length::new::<micrometer>(0.544)
}
/// returns thermal conductivity of copper
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
#[inline]
pub fn copper_specific_heat_capacity_zou_zweibaum_spline(
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
