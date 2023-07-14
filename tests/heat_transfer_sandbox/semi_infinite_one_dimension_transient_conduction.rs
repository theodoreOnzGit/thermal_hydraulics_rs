use std::f64::consts::PI;
use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::data_enum_structs::DataUserSpecifiedConvectionResistance;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::{link_heat_transfer_entity, HeatTransferInteractionType};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::SolidMaterial::SteelSS304L;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::Material;
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, OuterDiameterThermalConduction, SurfaceArea, SingleCVNode, CVType};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::CVType::SingleCV;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::density::density;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy};



use uom::si::f64::*;
use uom::si::length::centimeter;
use uom::si::power::watt;
use uom::si::temperature_interval::degree_celsius as interval_deg_c;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

/// This is an example of transient conduction case where analytical 
/// solutions have been well know 
///
/// This is conduction in a semi infinite medium with constant 
/// temperature boundary conditions 
///
/// These analytical solutions are well known
/// and contain error functions 
///
/// For constant temperature, the temperature at point x 
/// and time t is given by 
///
/// theta (x,t) = erfc (x / (2.0 * sqrt{alpha t}) )
///
///
/// Trojan, M. "Transient Heat Conduction in Semiinfinite 
/// Solid with Surface Convection." Encyclopedia of (2014).
///
/// One may note that the (x / (2.0 * sqrt{alpha t}) )
/// is simply the reciprocal of a Fourier like number multiplied by 
/// a constant
///
/// Fo = (alpha t)/L^2 
///
/// In this case, the length scale is x and t is the time.
///
/// How long must the medium be for it to be semi infinite? 
///
/// erfc (zeta) can be plotted on a graph 
/// at erfc(zeta) where zeta >= 2.0,
/// erfc(zeta) <= 0.005 (0.5\% change)
///
/// theta (x,t) = erfc (1 / (2.0 * sqrt{Fo(x,t)}) )
///
///
/// probably use some published reference, better quality
/// // https://www.unipamplona.edu.co/unipamplona/portalIG/home_34/recursos/01general/21082014/unidad_2_termo_ii.pdf
///
///
#[test]
fn transient_conduction_semi_infinite_copper_medium() -> Result<(), String>{

    // let's first do the analytical solution

    return Err("not finished".to_string());
}
