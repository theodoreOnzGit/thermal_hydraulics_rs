use uom::si::f64::*;

use crate::boussinesq_solver::pre_built_components::ciet_steady_state_natural_circulation_test_components::dracs_loop_components::*;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;

/// builds the hot branch of the DRACS loop (somewhat like the hot leg,
/// but with some other stuff)
/// this is pipe 30a all the way to 34 
/// but I build the components from top down
pub fn dracs_hot_branch_builder(initial_temperature: ThermodynamicTemperature) -> 
FluidComponentCollection {
    let pipe_34 = new_pipe_34(initial_temperature);
    let pipe_33 = new_pipe_33(initial_temperature);
    let pipe_32 = new_pipe_32(initial_temperature);
    let pipe_31a = new_pipe_31a(initial_temperature);
    let static_mixer_61_label_31 = new_static_mixer_61_label_31(initial_temperature);
    let dhx_tube_side_30b = new_dhx_tube_side_30b(initial_temperature);
    let dhx_tube_side_heat_exchanger_30 = new_isolated_dhx_tube_side_30(initial_temperature);
    let dhx_tube_side_30a = new_dhx_tube_side_30a(initial_temperature);



    let mut dracs_hot_branch = 
        FluidComponentCollection::new_series_component_collection();

    dracs_hot_branch.clone_and_add_component(&pipe_34);
    dracs_hot_branch.clone_and_add_component(&pipe_33);
    dracs_hot_branch.clone_and_add_component(&pipe_32);
    dracs_hot_branch.clone_and_add_component(&pipe_31a);
    dracs_hot_branch.clone_and_add_component(&static_mixer_61_label_31);
    dracs_hot_branch.clone_and_add_component(&dhx_tube_side_30b);
    dracs_hot_branch.clone_and_add_component(&dhx_tube_side_heat_exchanger_30);
    dracs_hot_branch.clone_and_add_component(&dhx_tube_side_30a);

    dracs_hot_branch
}


/// builds the cold branch of the DRACS loop (somewhat like the cold leg,
/// but with some other stuff)

pub fn dracs_cold_branch_builder(initial_temperature: ThermodynamicTemperature) -> 
FluidComponentCollection {
    let tchx_35a = new_inactive_ndhx_tchx_horizontal_35a(initial_temperature);
    let tchx_35b = new_ndhx_tchx_vertical_35b(initial_temperature);
    let static_mixer_60_label_36 = new_static_mixer_60_label_36(initial_temperature);
    let pipe_36a = new_pipe_36a(initial_temperature);
    let pipe_37 = new_pipe_37(initial_temperature);
    let flowmeter_60_37a = new_flowmeter_60_37a(initial_temperature);
    let pipe_38 = new_pipe_38(initial_temperature);
    let pipe_39 = new_pipe_39(initial_temperature);



    let mut dracs_cold_branch = 
        FluidComponentCollection::new_series_component_collection();

    dracs_cold_branch.clone_and_add_component(&tchx_35a);
    dracs_cold_branch.clone_and_add_component(&tchx_35b);
    dracs_cold_branch.clone_and_add_component(&static_mixer_60_label_36);
    dracs_cold_branch.clone_and_add_component(&pipe_36a);
    dracs_cold_branch.clone_and_add_component(&pipe_37);
    dracs_cold_branch.clone_and_add_component(&flowmeter_60_37a);
    dracs_cold_branch.clone_and_add_component(&pipe_38);
    dracs_cold_branch.clone_and_add_component(&pipe_39);

    dracs_cold_branch

}


#[test] 
pub fn dracs_branch_pressure_change_test(){

    // let's construct the branches with test pressures and obtain 
    // mass flowrates
    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_collection::FluidComponentCollectionMethods;
    use uom::si::f64::*;
    use uom::ConstZero;
    use approx::assert_abs_diff_eq;
    use uom::si::pressure::pascal;

    use uom::si::thermodynamic_temperature::degree_celsius;

    let test_temperature = ThermodynamicTemperature::
        new::<degree_celsius>(21.7);
    let mut dracs_hot_branch = dracs_hot_branch_builder(test_temperature);

    // pressure change at 0 kg/s 
    let pressure_change_at_zero_kg_per_s_hot_branch = 
        dracs_hot_branch.get_pressure_change(MassRate::ZERO);

    // pressure change should be 53627 +/- 1 Pa
    assert_abs_diff_eq!(pressure_change_at_zero_kg_per_s_hot_branch.get::<pascal>(), 
        53627.0, epsilon=1.0,);

    let dracs_cold_branch = dracs_cold_branch_builder(test_temperature);

    let pressure_change_at_zero_kg_per_s_cold_branch = 
        dracs_cold_branch.get_pressure_change(MassRate::ZERO);

    // pressure change should be same for both branches at equal flow
    assert_abs_diff_eq!(pressure_change_at_zero_kg_per_s_cold_branch.get::<pascal>(), 
        pressure_change_at_zero_kg_per_s_hot_branch.get::<pascal>(), 
        epsilon=1.0,);

    // now suppose I set the hot branch to 80 degrees C, there should 
    // be some buoyancy force thus, pressure change is less than 53627 Pa

    dracs_hot_branch = dracs_hot_branch_builder(
        ThermodynamicTemperature::new::<degree_celsius>(80.0));

    // pressure change at 0 kg/s 
    let pressure_change_at_zero_kg_per_s_hot_branch = 
        dracs_hot_branch.get_pressure_change(MassRate::ZERO);

    // pressure change should be less than 53627 +/- 1 Pa
    // in this case 51119 Pa
    assert_abs_diff_eq!(pressure_change_at_zero_kg_per_s_hot_branch.get::<pascal>(), 
        51119.0, epsilon=1.0,);
}
