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
