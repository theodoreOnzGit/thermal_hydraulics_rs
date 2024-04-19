use uom::{si::{mass_rate::kilogram_per_second, pressure::pascal}, ConstZero};

use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::{fluid_component_collection::FluidComponentCollectionMethods, fluid_component_super_collection::FluidComponentSuperCollection};


#[test]
pub fn heater_branch_with_heater_v2_test(){

    use crate::boussinesq_solver::pre_built_components::ciet_isothermal_test_components::{new_heated_section_version_1_label_1, new_heater_bottom_head_1b, new_heater_top_head_1a, new_pipe_18, new_pipe_2a, new_static_mixer_10};
    use uom::si::f64::*;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::{fluid_component::FluidComponent, fluid_component_collection::FluidComponentCollectionOreintation};
    use uom::si::mass_rate::kilogram_per_second;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollectionMethods;

    use uom::si::pressure::pascal;
    use super::{new_branch_5, new_pipe_3, new_pipe_4};
    // first let's construct the heater branch
    // probably need the heater top and bottom head later
    
    let branch_5 = new_branch_5();
    let pipe_4 = new_pipe_4();
    let pipe_3 = new_pipe_3();

    let static_mixer_2 = new_static_mixer_10();
    let static_mixer_pipe_2a = new_pipe_2a();
    let heater_top_head_1a = new_heater_top_head_1a();
    let heated_section_1 = new_heated_section_version_1_label_1();
    let heater_bottom_head_1b = new_heater_bottom_head_1b();

    let pipe_18 = new_pipe_18();

    // from top to bottom convention, that is branch 5 to 1b
    // but first, need to convert them into fluid components first 
    let branch_5_component: FluidComponent = 
        FluidComponent::FluidArray(
            branch_5.pipe_fluid_array.clone().try_into().unwrap()
            );

    let pipe_4_component: FluidComponent = 
        FluidComponent::FluidArray(
            pipe_4.pipe_fluid_array.clone().try_into().unwrap()
            );

    let pipe_3_component: FluidComponent = 
        FluidComponent::FluidArray(
            pipe_3.pipe_fluid_array.clone().try_into().unwrap()
            );

    let static_mixer_2_component: FluidComponent = 
        FluidComponent::FluidArray(
            static_mixer_2.pipe_fluid_array.try_into().unwrap()
            );

    let static_mixer_pipe_2a_component: FluidComponent = 
        FluidComponent::FluidArray(
            static_mixer_pipe_2a.pipe_fluid_array.try_into().unwrap()
            );

    let heater_top_head_1a_component: FluidComponent = 
        FluidComponent::FluidArray(
            heater_top_head_1a.pipe_fluid_array.try_into().unwrap()
            );

    let heater_1_component: FluidComponent = 
        FluidComponent::FluidArray(
            heated_section_1.pipe_fluid_array.try_into().unwrap()
            );

    let heater_bottom_head_1b_component: FluidComponent = 
        FluidComponent::FluidArray(
            heater_bottom_head_1b.pipe_fluid_array.try_into().unwrap()
            );

    let pipe_18_component: FluidComponent = 
        FluidComponent::FluidArray(
            pipe_18.pipe_fluid_array.try_into().unwrap()
            );


    let heater_branch = 
        FluidComponentCollection{
            components: vec![
                branch_5_component,
                pipe_4_component,
                pipe_3_component,
                static_mixer_2_component,
                static_mixer_pipe_2a_component,
                heater_top_head_1a_component,
                heater_1_component,
                heater_bottom_head_1b_component,
                pipe_18_component
            ],
            orientation: FluidComponentCollectionOreintation::Series,
        };


    // let's check the hydrostatic pressure, 0.0 kg/s fluid flow 
    {
        // now let's push a 0.1kg/s fluid flow through this pipe series
        //
        let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.0);

        // and then let's get the pressure change

        let series_pipe_pressure_change = heater_branch.
            get_pressure_change(pipe_fluid_flow);

        // pressure change is around 39041 Pa
        approx::assert_relative_eq!(
            series_pipe_pressure_change.get::<pascal>(),
            39041.0,
            max_relative=0.001);
    }
    // now let's push a 0.1kg/s fluid flow through this pipe series
    //
    let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.1);

    // and then let's get the pressure change

    let series_pipe_pressure_change = heater_branch.
        get_pressure_change(pipe_fluid_flow);

    // pressure change is around 36451 Pa
    approx::assert_relative_eq!(
        series_pipe_pressure_change.get::<pascal>(),
        36451.0,
        max_relative=0.001);
}


#[test]
pub fn ctah_branch_test(){

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

    // let's check the hydrostatic pressure, 0.0 kg/s fluid flow 
    {
        // now let's push a 0.1kg/s fluid flow through this pipe series
        //
        let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.0);

        // and then let's get the pressure change

        let series_pipe_pressure_change = ctah_branch.
            get_pressure_change(pipe_fluid_flow);

        // pressure change is around 39041 Pa
        approx::assert_relative_eq!(
            series_pipe_pressure_change.get::<pascal>(),
            39041.0,
            max_relative=0.001);
    }
}


#[test]
pub fn isothermal_ctah_and_heater_branch_validation_test(){

    // let's construct the branches with test pressures and obtain 
    // mass flowrates
    use uom::si::f64::*;
    use super::ciet_branch_builders_isothermal::{ctah_branch_builder_isothermal_test, heater_branch_builder_isothermal_test};

    {
        let test_pressure_1 = Pressure::new::<pascal>(0.0);
        let heater_branch = heater_branch_builder_isothermal_test();
        let ctah_branch = ctah_branch_builder_isothermal_test(test_pressure_1);

        // expected flowrate 
        let expected_mass_flow_1 = MassRate::new::<kilogram_per_second>(0.0);

        // you'll now need to add this into a super collection
        let mut ctah_and_heater_branches = 
            FluidComponentSuperCollection::default();
        ctah_and_heater_branches.set_oritentation_to_parallel();
        ctah_and_heater_branches.fluid_component_super_vector
            .push(heater_branch);
        ctah_and_heater_branches.fluid_component_super_vector
            .push(ctah_branch);

        // obtain flowrate 
        // the overall pressure change across each branch must be
        // equal
        //
        // the total mass flowrate through the collection of parallel 
        // branches is zero
        let net_mass_flowrate_through_parallel_branches = 
            MassRate::ZERO;
        let pressure_change_across_each_branch = 
            ctah_and_heater_branches.get_pressure_change(
            net_mass_flowrate_through_parallel_branches);

        // once we have pressure change across each branch,
        // then we can calculate mass flowrate.

        let mass_flowrate_across_each_branch: Vec<MassRate> = 
            ctah_and_heater_branches.
            get_mass_flowrate_across_each_parallel_branch(
                pressure_change_across_each_branch
            );




    }


}


