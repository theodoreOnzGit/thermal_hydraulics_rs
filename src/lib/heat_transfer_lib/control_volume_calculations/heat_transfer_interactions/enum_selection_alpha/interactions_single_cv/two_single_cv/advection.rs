
use enums_alpha::data_enum_structs::DataAdvection;
use uom::num_traits::Zero;
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
    advection_data: DataAdvection)-> Result<(), String>{

    let mass_flow_from_cv_1_to_cv_2 = advection_data.mass_flowrate;

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
    //
    // the timescale can only be calculated after the mass flows 
    // in and out of the cv are sufficiently calculated
    // the only thing we can do here is push the mass flowrate 
    // into the individual mass flowrate vectors 
    //
    // by convention, mass flowrate goes out of cv1 and into cv2 
    // so positive mass flowrate here is positive for cv2 
    // and negative for cv1
    //
    // I'll need a density for the flow first

    let density_cv1 = advection_data.fluid_density_cv1;
    let density_cv2 = advection_data.fluid_density_cv2;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_cv_1_to_cv_2 > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from cv1 
        // to cv2 
        // then the density we use is cv1 

        volumetric_flowrate = mass_flow_from_cv_1_to_cv_2/density_cv1;

    } else {
        // if mass flowrate is positive, flow is moving from cv2
        // to cv1
        // then the density we use is cv2

        volumetric_flowrate = mass_flow_from_cv_1_to_cv_2/density_cv2;
    }

    // now that I've done the volume flowrate calculation, push the 
    // volumetric flowrate to each vector
    single_cv_1.volumetric_flowrate_vector.push(
        -volumetric_flowrate);
    single_cv_2.volumetric_flowrate_vector.push(
        volumetric_flowrate);

    

    // done! 
    Ok(())
}
