/// handles cases where two cvs have heat transfer governed by 
/// thermal conductance or thermal resistance 
pub mod conductance_interactions;
pub (in crate) use conductance_interactions::*;

/// handles cases where two cvs have heat transfer governed by 
/// advection or fluid flow 
pub mod advection;
pub (in crate) use advection::*;

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;

/// this is mostly a wrapper function
///
/// which calls other functions depending on whether the 
/// heat transfer interaction is conductance based on advection based
#[inline]
pub fn calculate_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType)-> Result<(), String>{


    match interaction {
        HeatTransferInteractionType::UserSpecifiedThermalConductance(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::SingleCartesianThermalConductanceOneDimension(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCartesianThermalConductanceThreeDimension(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCartesianThermalConductance(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCylindricalThermalConductance(_, _, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::CylindricalConductionConvectionLiquidOutside(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::CylindricalConductionConvectionLiquidInside(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatAddition => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCustomArea(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCylindricalOuterArea(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCylindricalInnerArea(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedConvectionResistance(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
    }

}

