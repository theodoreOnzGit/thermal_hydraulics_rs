use uom::si::f64::*;
use uom::si::length::meter;

use crate::boussinesq_solver::pre_built_components::shell_and_tube_heat_exchanger::SimpleShellAndTubeHeatExchanger;
use crate::prelude::beta_testing::SolidMaterial;

/// constructs a new instance of the shell and tube 
/// heat exchanger for the DHX based on Zou's specifications
/// for the flow area, hydraulic diameter and number of tubes
/// but Zweibaum's specifications for insulation thickness 
/// 
///
/// the heat transfer coefficients are based on Gnielinski 
/// correlation 
///
/// Whereas hydrodynamically, the DHX shell and tube 
/// sides are modelled as pipes with K values of 23.9 on 
/// the tube side and 3.3 on the tube side 
/// insulation thickness for DHX is 0.0508 m of fiberglass
/// DHX is made of copper tubing on the inside
/// and assumed to be copper on shell side as well
pub fn new_dhx_sthe_version_1() -> SimpleShellAndTubeHeatExchanger {

    let insulation_thickness: Length = Length::new::<meter>(0.0508);
    let copper = SolidMaterial::Copper;

    todo!()
}
