
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;



use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
#[inline]
pub fn calculate_advection_interaction_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    mass_flow_from_cv_1_to_cv_2: MassRate)-> Result<(), String>{

    // for this, quite straightforward, 
    // get both specific enthalpy of both cvs 
    // calculate power flow 
    // and then update the power vector 
    //

    let specific_enthalpy_cv1: AvailableEnergy = 
    single_cv_1.current_timestep_control_volume_specific_enthalpy;

    let specific_enthalpy_cv2: AvailableEnergy = 
    single_cv_2.current_timestep_control_volume_specific_enthalpy;

    // calculate heat rate 

    let heat_flowrate_from_cv_1_to_cv_2: Power 
    = advection_heat_rate(mass_flow_from_cv_1_to_cv_2,
        specific_enthalpy_cv1,
        specific_enthalpy_cv2,)?;

    // by default, cv 1 is on the left, cv2 is on the right 
    //

    single_cv_1.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_cv_1_to_cv_2);
    single_cv_2.rate_enthalpy_change_vector.
        push(heat_flowrate_from_cv_1_to_cv_2);

    // relevant timescale here is courant number

    // done! 
    Ok(())
}
