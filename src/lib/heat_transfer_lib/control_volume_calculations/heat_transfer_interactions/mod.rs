//!  For heat_transfer_interactions, 
//!
//!  the way this library is intended to be used is that 
//!  the user creates two heat transfer entities, they could 
//!  be control volumes or boundary conditions, 
//!  and then connects them using some thermal resistance 
//!  or thermal conductance
//!
//! The calculation logic is stored in a calculation file 
//!
//! But the selection logic is strongly based on enums, results and 
//! typing hopefully so that the calculations would be more or less 
//! correct at compiletime rather than runtime
//!
//!




use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::HeatTransferEntity;

/// This is a calculation library for all interactions between 
/// the various heat transfer entities
///
/// but does not contain the selection logic yet
///
/// Stuff here is not included in the API, so it might change
pub(in crate::heat_transfer_lib:: control_volume_calculations)
mod calculations;
use calculations::*;


/// contains enums for the user to select how heat transfer entities 
/// ie control volumes and boundary conditions 
/// interact with each other
///
/// note that the enums here are UNSTABLE, I may not preserve the API
/// alpha denotes unstable, like program alpha and beta testing
pub mod enums_alpha;
pub use enums_alpha::*;

/// contains scripts to select enums 
/// for how heat transfer entities 
/// ie control volumes and boundary conditions 
/// interact with each other
/// and also control volume mutation
/// 
///
/// one thread calculation only
///
/// note that the enum_selection_alpha module 
/// here are UNSTABLE, I may not preserve the API
///
/// alpha here means unstable 
///
///
pub mod enum_selection_alpha;
use enum_selection_alpha::*;
use uom::si::f64::*;



///// contains scripts to select enums 
///// for how heat transfer entities 
///// ie control volumes and boundary conditions 
///// interact with each other
///// and also control volume mutation
/////
///// multithread calculation only
/////
///// note that the enum_selection_alpha module 
///// here are UNSTABLE, I may not preserve the API
/////
///// alpha here means unstable
////mod enum_selection_parallel_alpha;
////use enum_selection_parallel_alpha::*;


/// For this part, we determine how heat_transfer_entities interact 
/// with each other by using a function 
/// The function will take in a HeatTransferInteractionType enum
/// which you must first initiate
///
/// Then you need to supply two control volumes or more generally 
/// heat_transfer_entities, which can consist of mix of control volumes
/// and boundary conditions
///
/// The function will then calculate the heat transfer between the two 
/// control volumes, and either return a value or mutate the CV objects 
/// using mutable borrows
pub fn link_heat_transfer_entity(entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType)-> Result<(), String>{

    // first thing first, probably want to unpack the enums to obtain 
    // the underlying control volume and BCs
    // 
    // Basically there are four permutations of what the user may choose 
    // either link two control volumes,
    // link a control vol and BC  (or BC first then control vol)
    // lastly, two BCs 
    // (which is kind of invalid though, 
    // but i suppose it may make sense for steady state)
    //

    let heat_transfer_entity_result = match (entity_1,entity_2) {
        (HeatTransferEntity::ControlVolume(cv_type_1),
            HeatTransferEntity::ControlVolume(cv_type_2)) =>
            calculate_control_volume_serial(
                cv_type_1, 
                cv_type_2, 
                interaction),
        (HeatTransferEntity::BoundaryConditions(bc_type),
            HeatTransferEntity::ControlVolume(cv_type)) =>
            attach_boundary_condition_to_control_volume_back_serial(
                bc_type, cv_type, interaction),
        (HeatTransferEntity::ControlVolume(cv_type),
            HeatTransferEntity::BoundaryConditions(bc_type)) =>
            attach_boundary_condition_to_control_volume_front_serial(
                cv_type, bc_type, interaction),
        (HeatTransferEntity::BoundaryConditions(bc_type_1),
            HeatTransferEntity::BoundaryConditions(bc_type_2)) =>
            calculate_boundary_condition_serial(
                bc_type_1, bc_type_2, interaction),
    };


    return heat_transfer_entity_result;
   
}



/// this function calculates relevant timescales when linking 
/// two heat transfer entities
pub fn calculate_timescales_for_heat_transfer_entity(
    entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType)-> Result<Time, String>{

    // first thing first, probably want to unpack the enums to obtain 
    // the underlying control volume and BCs
    // 
    // Basically there are four permutations of what the user may choose 
    // either link two control volumes,
    // link a control vol and BC  (or BC first then control vol)
    // lastly, two BCs 
    // (which is kind of invalid though, 
    // but i suppose it may make sense for steady state)
    //

    let heat_transfer_entity_timestep_result = match (entity_1,entity_2) {
        (HeatTransferEntity::ControlVolume(cv_type_1),
            HeatTransferEntity::ControlVolume(cv_type_2)) =>
            {
                calculate_timestep_control_volume_serial(
                            cv_type_1, 
                            cv_type_2, 
                            interaction)
            },
        (HeatTransferEntity::BoundaryConditions(bc_type),
            HeatTransferEntity::ControlVolume(cv_type)) =>
            {
                calculate_timestep_control_volume_back_to_boundary_condition_serial(
                    bc_type, cv_type, interaction)
            },
        (HeatTransferEntity::ControlVolume(cv_type),
            HeatTransferEntity::BoundaryConditions(bc_type)) =>
            {
                calculate_timestep_control_volume_front_to_boundary_condition_serial(
                    cv_type, bc_type, interaction)
            },
        (HeatTransferEntity::BoundaryConditions(bc_type_1),
            HeatTransferEntity::BoundaryConditions(bc_type_2)) =>
            {
                calculate_timestep_boundary_condition_serial(
                    bc_type_1,
                    bc_type_2,
                    interaction)
            }
    };


    return heat_transfer_entity_timestep_result;
   
}



