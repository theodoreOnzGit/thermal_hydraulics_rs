
use uom::si::f64::*;
use crate::heat_transfer_lib;



use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
#[inline]
pub fn calculate_single_cv_node_constant_heat_addition(
    control_vol: &mut SingleCVNode,
    heat_added_to_control_vol: Power,
    interaction: HeatTransferInteractionType
    ) -> Result<(), String> {

    // ensure that the interaction is UserSpecifiedHeatAddition
    // otherwise, return error 

    match interaction {
        heat_transfer_lib::control_volume_calculations::
            heat_transfer_interactions::
            enums_alpha::HeatTransferInteractionType::
            UserSpecifiedHeatAddition => {
                // return a void value, that would be dropped 
                // instantly
                //
                // it pretty much has the same meaning as break
                ()
            },

        _ => return Err("you need to specify that the interaction type \n 
            is UserSpecifiedHeatAddition".to_string()),
    };


    control_vol.rate_enthalpy_change_vector.
        push(heat_added_to_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change
    //
    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    


    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;

            ()
        },
        Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }

    return Ok(());
}
