





/// Here is a test which is meant to test a simple struct made
/// to hold and calculate fluid component collections
///
/// First i make a typical fluid component, a set of air pipes
/// perhaps 10 air pipes and i want to put them in series
#[test]
pub fn simple_fluid_collection_example_5 () {

    // this tests the calc pressure loss for fluid component 
    use uom::si::f64::*;
    use uom::ConstZero;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidCustomComponentCalcPressureLoss;
    use uom::si::length::inch;
    use std::f64::consts::PI;
    use uom::si::ratio::ratio;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::mass_density::kilogram_per_cubic_meter;
    use uom::si::dynamic_viscosity::poise;
    use uom::si::pressure::pascal;
    use uom::si::length::meter;
    use uom::si::angle::degree;
    use uom::si::length::millimeter;
    use uom::si::dynamic_viscosity::millipascal_second;

    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidComponentTrait;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidCustomComponentCalcPressureChange;
    use crate::boussinesq_solver::fluid_mechanics_correlations::pipe_calculations::pipe_calc_mass_flowrate;
    use crate::boussinesq_solver::fluid_mechanics_correlations::pipe_calculations::pipe_calc_pressure_loss;

    use uom::si::{pressure::atmosphere, thermodynamic_temperature::kelvin};

    use crate::boussinesq_solver::boussinesq_thermophysical_properties::{LiquidMaterial, SolidMaterial};
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component::FluidComponent;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollectionOreintation;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_collection::FluidComponentCollection;


    // honestly, I am quite unsure as to how to structure the 
    //
    // fluid mechanics and thermal hydraulics calculations.
    // 
    //
    // there are several ways we can do this... 
    //
    //
    // We could first manage one component at a time, and manually 
    // sum up the pressure drop contributions of each component as well 
    // as manually update each of the component temperatures
    //
    // This, however, is a very cumbersome way to coding. It is extremely 
    // error prone.
    //
    // To abstract away some of this, we could have individual heat 
    // transfer entities, and at each timestep, clone their relevant fluid 
    // arrays into a FluidComponentCollection. These FluidComponentCollections 
    // can then be arranged into a FluidComponentSuperCollection
    // and then the relevant fluid mechanics calculations could be 
    // performed.
    //
    // The heat transfer entities would then have their 
    // mass flowrates sorted out manually at each timestep.
    //
    // The most abstract way would be to couple all components into a 
    // collection. However, given that each component has multiple 
    // heat transfer entities in it, this may not be the best way to do 
    // things.
    //
    // For now, let's opt for way two, where we have HeatTransferEntity 
    // objects converted into FluidComponent objects, and then arranged 
    // into arrays.
    //
    // HeatTransferEntity objects are usually in the pre_built_components 
    // section. So I won't test them here.
    // However, these can be converted into FluidComponent objects.
    // So we'll just work with fluid component objects for now.

    let adjacent_solid_material = SolidMaterial::Copper;
    let liquid_material = LiquidMaterial::TherminolVP1;
    let initial_pressure = Pressure::new::<atmosphere>(1.0);
    let initial_temperature = ThermodynamicTemperature::new::<kelvin>(298.0);
    let hydraulic_diameter = Length::new::<inch>(2.0);
    let pipe_incline_angle = Angle::new::<degree>(0.0);
    let pipe_form_loss = 5.0;
    let user_specified_inner_nodes = 0;
    let length = Length::new::<meter>(1.0);

    let therminol_pipe_1: FluidComponent = FluidComponent::FluidArray(
        FluidArray::new_cylinder(
            length, 
            hydraulic_diameter, 
            initial_temperature, 
            initial_pressure, 
            adjacent_solid_material, 
            liquid_material, 
            pipe_form_loss.into(), 
            user_specified_inner_nodes, 
            pipe_incline_angle)
        );

    let therminol_pipe_2 = therminol_pipe_1.clone();
    let therminol_pipe_3 = therminol_pipe_1.clone();
    let therminol_pipe_4 = therminol_pipe_1.clone();
    let therminol_pipe_5 = therminol_pipe_1.clone();
    let therminol_pipe_6 = therminol_pipe_1.clone();
    let therminol_pipe_7 = therminol_pipe_1.clone();
    let therminol_pipe_8 = therminol_pipe_1.clone();
    let therminol_pipe_9 = therminol_pipe_1.clone();
    let therminol_pipe_10 = therminol_pipe_1.clone();

    // let's now put all of them in series.
    
    let series_collection_of_therminol_pipes = 
        FluidComponentCollection{
            components: vec![therminol_pipe_1.clone(),
            therminol_pipe_2.clone(),
            therminol_pipe_3.clone(),
            therminol_pipe_4.clone(),
            therminol_pipe_5.clone(),
            therminol_pipe_6.clone(),
            therminol_pipe_7.clone(),
            therminol_pipe_8.clone(),
            therminol_pipe_9.clone(),
            therminol_pipe_10.clone(),],
            orientation: FluidComponentCollectionOreintation::Series,
        };






}
