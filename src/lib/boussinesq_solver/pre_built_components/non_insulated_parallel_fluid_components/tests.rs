#[test] 
pub fn parallel_tube_mass_flowrate_should_be_more_than_single_tube(){

    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30_parallel_tubes;
    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30;
    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::pascal;

    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_traits::FluidComponentTrait;

    let initial_temperature = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);

    let parallel_dhx = 
        new_isolated_dhx_tube_side_30_parallel_tubes(initial_temperature);
    let single_dhx = 
        new_isolated_dhx_tube_side_30(initial_temperature);

    let pressure_loss = Pressure::new::<pascal>(1000.0);

    // mass flowrates 
    let single_dhx_massrate: MassRate = 
        single_dhx.get_mass_flowrate_from_pressure_loss_immutable(pressure_loss);

    approx::assert_abs_diff_eq!(
        single_dhx_massrate.get::<kilogram_per_second>(),
        0.275,
        epsilon=0.001
        );

    let parallel_dhx_mass_flowrate: MassRate = 
        parallel_dhx.get_mass_flowrate_from_pressure_loss_immutable(pressure_loss);

    approx::assert_abs_diff_eq!(
        parallel_dhx_mass_flowrate.get::<kilogram_per_second>(),
        0.275 * (parallel_dhx.number_of_tubes as f64),
        epsilon=0.01
        );
}


#[test] 
pub fn parallel_tube_set_mass_flowrate_overall_should_correspond_to_single_tube(){

    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30_parallel_tubes;
    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::mass_rate::kilogram_per_second;

    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_traits::FluidComponentTrait;
    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        one_d_fluid_array_with_lateral_coupling::FluidArray;

    let initial_temperature = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);

    let mut parallel_dhx = 
        new_isolated_dhx_tube_side_30_parallel_tubes(initial_temperature);

    // we have 19 tubes 
    // so the total mass flowrate across 19 tubes is set at 0.19 kg/s 
    //
    //
    let total_mass_flowrate = MassRate::new::<kilogram_per_second>(0.19);

    // set and get mass flowrates 

    parallel_dhx.set_mass_flowrate(total_mass_flowrate);


    approx::assert_abs_diff_eq!(
        parallel_dhx.get_mass_flowrate().get::<kilogram_per_second>(),
        0.19,
        epsilon=0.001
        );

    // now let's get the mass flowrate for the single tube 

    let mut pipe_fluid_array: FluidArray = 
        parallel_dhx.pipe_fluid_array.clone().try_into().unwrap();


    // the mass flowrate should be 0.01 kg/s
    approx::assert_abs_diff_eq!(
        pipe_fluid_array.get_mass_flowrate().get::<kilogram_per_second>(),
        0.01,
        epsilon=0.001
        );


}


#[test] 
pub fn parallel_tube_get_pressure_loss_from_mass_flowrate(){

    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30_parallel_tubes;
    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30;
    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::pascal;

    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_traits::FluidComponentTrait;

    let initial_temperature = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);

    let parallel_dhx = 
        new_isolated_dhx_tube_side_30_parallel_tubes(initial_temperature);

    let reference_pressure_loss = Pressure::new::<pascal>(1000.0);
    // for this, we have 19 tubes
    let reference_mass_flowrate = MassRate::new::<kilogram_per_second>(0.275*19.0);

    let parallel_dhx_pressure_loss: Pressure = 
        parallel_dhx.get_pressure_loss_immutable(reference_mass_flowrate);

    approx::assert_abs_diff_eq!(
        parallel_dhx_pressure_loss.get::<pascal>(),
        reference_pressure_loss.get::<pascal>(),
        epsilon=5.0 // max difference of 5 Pa
        );
}
