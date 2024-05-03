

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

    let mut parallel_dhx = 
        new_isolated_dhx_tube_side_30_parallel_tubes(initial_temperature);
    let mut single_dhx = 
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
