use uom::si::pressure::pascal;








#[test]
pub fn heater_branch_with_heater_v2_test(){

    use crate::boussinesq_solver::pre_built_components::ciet_static_mixers::StaticMixers;
    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use crate::boussinesq_solver::pre_built_components::ciet_heater_top_and_bottom_head_bare::HeaterTopBottomHead;
    use crate::boussinesq_solver::pre_built_components::ciet_heater_version_2_bare::HeaterVersion2Bare;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::{fluid_component::FluidComponent, fluid_component_collection::FluidComponentCollectionOreintation};
    use uom::si::mass_rate::kilogram_per_second;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollectionMethods;

    use super::{new_branch_5, new_pipe_3, new_pipe_4};
    // first let's construct the heater branch
    // probably need the heater top and bottom head later
    
    let branch_5 = new_branch_5();
    let pipe_4 = new_pipe_4();
    let pipe_3 = new_pipe_3();

    let initial_temperature = ThermodynamicTemperature::new::<
        degree_celsius>(21.7);
    let ambient_temperature = ThermodynamicTemperature::new::<
        degree_celsius>(20.0);
    let static_mixer_2 = StaticMixers::new_static_mixer_2_mx10(
        initial_temperature, ambient_temperature);
    let static_mixer_pipe_2a = StaticMixers::new_static_mixer_pipe_2a_mx10(
        initial_temperature, ambient_temperature);
    // placeholders for now

    let heater_top_head_1a = HeaterTopBottomHead::new_top_head(
        initial_temperature, ambient_temperature);
    let heated_section_1 = HeaterVersion2Bare::new_dewet_model(
        initial_temperature, ambient_temperature, 6);
    let heater_bottom_head_1b = HeaterTopBottomHead::new_bottom_head(
        initial_temperature, ambient_temperature);

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
            static_mixer_2.therminol_array.try_into().unwrap()
            );

    let static_mixer_pipe_2a_component: FluidComponent = 
        FluidComponent::FluidArray(
            static_mixer_pipe_2a.therminol_array.try_into().unwrap()
            );

    let heater_top_head_1a_component: FluidComponent = 
        FluidComponent::FluidArray(
            heater_top_head_1a.therminol_array.try_into().unwrap()
            );

    let heater_1_component: FluidComponent = 
        FluidComponent::FluidArray(
            heated_section_1.therminol_array.try_into().unwrap()
            );

    let heater_bottom_head_1b_component: FluidComponent = 
        FluidComponent::FluidArray(
            heater_bottom_head_1b.therminol_array.try_into().unwrap()
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
                heater_bottom_head_1b_component
            ],
            orientation: FluidComponentCollectionOreintation::Series,
        };

    // now let's push a 0.1kg/s fluid flow through this pipe series
    //
    let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.1);

    // and then let's get the pressure change

    let series_pipe_pressure_change = heater_branch.
        get_pressure_change(pipe_fluid_flow);

    // pressure change is around -10283 Pa
    approx::assert_relative_eq!(
        series_pipe_pressure_change.get::<pascal>(),
        -10283.0,
        max_relative=0.001);

    // let's check the hydrostatic pressure, 0.0 kg/s fluid flow 
    {
        // now let's push a 0.1kg/s fluid flow through this pipe series
        //
        let pipe_fluid_flow = MassRate::new::<kilogram_per_second>(0.0);

        // and then let's get the pressure change

        let series_pipe_pressure_change = heater_branch.
            get_pressure_change(pipe_fluid_flow);

        // pressure change is around -9735 Pa
        approx::assert_relative_eq!(
            series_pipe_pressure_change.get::<pascal>(),
            -9735.0,
            max_relative=0.001);
    }
}
