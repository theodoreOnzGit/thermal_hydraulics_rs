use uom::si::mass_rate::kilogram_per_second;

use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::{fluid_component_collection::FluidComponentCollection, fluid_component_traits::FluidComponentTrait};

use super::*;


/// builds a heater branch to simulate isothermal testing of ciet
pub fn heater_branch_builder_isothermal_test() -> FluidComponentCollection {

    let branch_5 = new_branch_5();
    let pipe_4 = new_pipe_4();
    let pipe_3 = new_pipe_3();
    let static_mixer_2 = new_static_mixer_10();
    let static_mixer_pipe_2a = new_pipe_2a();
    let heater_top_head_1a = new_heater_top_head_1a();
    let heated_section_1 = new_heated_section_version_1_label_1();
    let heater_bottom_head_1b = new_heater_bottom_head_1b();
    let pipe_18 = new_pipe_18();

    let mut heater_branch = FluidComponentCollection::new_series_component_collection();

    heater_branch.clone_and_add_component(&branch_5);
    heater_branch.clone_and_add_component(&pipe_4);
    heater_branch.clone_and_add_component(&pipe_3);
    heater_branch.clone_and_add_component(&static_mixer_2);
    heater_branch.clone_and_add_component(&static_mixer_pipe_2a);
    heater_branch.clone_and_add_component(&heater_top_head_1a);
    heater_branch.clone_and_add_component(&heated_section_1);
    heater_branch.clone_and_add_component(&heater_bottom_head_1b);
    heater_branch.clone_and_add_component(&pipe_18);

    heater_branch
}

/// builds the ctah branch to simulate isothermal testing of ciet 
/// allows user to supply a pump pressure or loop pressure drop
pub fn ctah_branch_builder_isothermal_test(
    pump_pressure: Pressure) -> FluidComponentCollection {

    let static_mixer_41_label_6 = new_static_mixer_41();
    let pipe_6a = new_pipe_6a();
    let ctah_vertical_label_7a = new_inactive_ctah_vertical();
    let ctah_horizontal_label_7b = new_inactive_ctah_horizontal();
    let pipe_8a = new_pipe_8a();
    let static_mixer_40_label_8 = new_static_mixer_40();
    let pipe_9 = new_pipe_9();
    let pipe_10 = new_pipe_10();
    let pipe_11 = new_pipe_11();
    let pipe_12 = new_pipe_12();
    let mut ctah_pump = new_ctah_pump();
    ctah_pump.set_internal_pressure_source(pump_pressure);
    let pipe_13 = new_pipe_13();
    let pipe_14 = new_pipe_14();
    let flowmeter_40_14a = new_flowmeter_40_14a();
    let pipe_15 = new_pipe_15();
    let pipe_16 = new_pipe_16();
    let branch_17 = new_branch_17();


    // now I want to add each of these to the fluid component 
    // collection without the constant hassle of having to convert types

    let mut ctah_branch = FluidComponentCollection::new_series_component_collection();
    
    ctah_branch.clone_and_add_component(&static_mixer_41_label_6);
    ctah_branch.clone_and_add_component(&pipe_6a);
    ctah_branch.clone_and_add_component(&ctah_vertical_label_7a);
    ctah_branch.clone_and_add_component(&ctah_horizontal_label_7b);
    ctah_branch.clone_and_add_component(&pipe_8a);
    ctah_branch.clone_and_add_component(&static_mixer_40_label_8);
    ctah_branch.clone_and_add_component(&pipe_9);
    ctah_branch.clone_and_add_component(&pipe_10);
    ctah_branch.clone_and_add_component(&pipe_11);
    ctah_branch.clone_and_add_component(&pipe_12);
    ctah_branch.clone_and_add_component(&ctah_pump);
    ctah_branch.clone_and_add_component(&pipe_13);
    ctah_branch.clone_and_add_component(&pipe_14);
    ctah_branch.clone_and_add_component(&flowmeter_40_14a);
    ctah_branch.clone_and_add_component(&pipe_15);
    ctah_branch.clone_and_add_component(&pipe_16);
    ctah_branch.clone_and_add_component(&branch_17);

    ctah_branch
}

/// for the CTAH pump, I expect zero resistance or pressure drop 
/// that is assumed
///
/// This is correct
#[test]
pub fn ctah_pump_should_give_zero_resistance(){
    let ctah_pump = new_ctah_pump();

    let mass_rate = MassRate::new::<kilogram_per_second>(0.18);

    let pressure_drop = ctah_pump.get_pressure_loss_immutable(mass_rate);

    approx::assert_abs_diff_eq!(
        pressure_drop.get::<uom::si::pressure::pascal>(),
        0.0,
        );

    let pressure_change = ctah_pump.get_pressure_change_immutable(mass_rate);

    approx::assert_abs_diff_eq!(
        pressure_change.get::<uom::si::pressure::pascal>(),
        0.0,
        )

}
