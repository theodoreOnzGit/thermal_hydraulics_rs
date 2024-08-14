use uom::si::f64::*;
use uom::si::length::millimeter;
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
pub fn fiberglass_density() -> Result<MassDensity,ThermalHydraulicsLibError> {
    return Ok(MassDensity::new::<kilogram_per_cubic_meter>(20.0));
}

/// Value from: Perry's chemical Engineering handbook 
/// 8th edition Table 6-1 
/// generic value for drawn tubing
/// Perry, R. H., & DW, G. (2007). 
/// Perry’s chemical engineers’ handbook, 
/// 8th illustrated ed. New York: McGraw-Hill.
pub fn fiberglass_surf_roughness() -> Length {
    Length::new::<millimeter>(0.00152)
}
