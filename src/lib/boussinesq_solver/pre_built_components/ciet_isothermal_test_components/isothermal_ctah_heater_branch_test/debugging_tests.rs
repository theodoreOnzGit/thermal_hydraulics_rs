/// basically the ctah branch is rather buggy, so we have to 
/// test one branch at a time
#[test]
pub fn partial_ctah_branch_test(){

    use crate::boussinesq_solver::pre_built_components::ciet_isothermal_test_components::{new_branch_17, new_ctah_pump, new_flowmeter_40_14a, new_inactive_ctah_horizontal, new_inactive_ctah_vertical, new_pipe_10, new_pipe_11, new_pipe_12, new_pipe_13, new_pipe_14, new_pipe_15, new_pipe_16, new_pipe_6a, new_pipe_8a, new_pipe_9, new_static_mixer_40, new_static_mixer_41};
    use uom::si::f64::*;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;
    use uom::si::mass_rate::kilogram_per_second;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollectionMethods;

    use uom::si::pressure::pascal;
    // first let's construct the ctah branch
    // this is pipe 6 all the way to branch 17

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
    let ctah_pump = new_ctah_pump();
    let pipe_13 = new_pipe_13();
    let pipe_14 = new_pipe_14();
    let flowmeter_40_14a = new_flowmeter_40_14a();
    let pipe_15 = new_pipe_15();
    let pipe_16 = new_pipe_16();
    let branch_17 = new_branch_17();


    // now I want to add each of these to the fluid component 
    // collection without the constant hassle of having to convert types
    //
    // first let's test up to pipe 10

    let mut ctah_branch = FluidComponentCollection::new_series_component_collection();
    
    ctah_branch.clone_and_add_component(&static_mixer_41_label_6);
    ctah_branch.clone_and_add_component(&pipe_6a);
    ctah_branch.clone_and_add_component(&ctah_vertical_label_7a);
    ctah_branch.clone_and_add_component(&ctah_horizontal_label_7b);
    ctah_branch.clone_and_add_component(&pipe_8a);
    ctah_branch.clone_and_add_component(&static_mixer_40_label_8);
    ctah_branch.clone_and_add_component(&pipe_9);
    ctah_branch.clone_and_add_component(&pipe_10);
    //ctah_branch.clone_and_add_component(&pipe_11);
    //ctah_branch.clone_and_add_component(&pipe_12);
    //ctah_branch.clone_and_add_component(&ctah_pump);
    //ctah_branch.clone_and_add_component(&pipe_13);
    //ctah_branch.clone_and_add_component(&pipe_14);
    //ctah_branch.clone_and_add_component(&flowmeter_40_14a);
    //ctah_branch.clone_and_add_component(&pipe_15);
    //ctah_branch.clone_and_add_component(&pipe_16);
    //ctah_branch.clone_and_add_component(&branch_17);

    // let's the pressure change, 0.18 kg/s fluid flow 
    {
        // now let's push a 0.18kg/s fluid flow through this pipe series
        //
        let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.18);

        // and then let's get the pressure change

        let series_pipe_pressure_change = ctah_branch.
            get_pressure_change(pipe_fluid_flow);

        // pressure change is around 39041 Pa
        approx::assert_relative_eq!(
            series_pipe_pressure_change.get::<pascal>(),
            28347.0,
            max_relative=0.001);
    }
}
