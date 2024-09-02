use uom::si::f64::*;

use crate::boussinesq_solver::
array_control_vol_and_fluid_component_collections::
fluid_component_collection::
fluid_component_collection::FluidComponentCollection;
// let's construct the branches with test pressures and obtain 
use crate::boussinesq_solver::
array_control_vol_and_fluid_component_collections::
fluid_component_collection::
fluid_component_collection::FluidComponentCollectionMethods;
use uom::ConstZero;

use uom::si::thermodynamic_temperature::degree_celsius;
use crate::boussinesq_solver::
array_control_vol_and_fluid_component_collections::
fluid_component_collection::
fluid_component_super_collection::FluidComponentSuperCollection;

use crate::boussinesq_solver::pre_built_components::
insulated_pipes_and_fluid_components::InsulatedFluidComponent;
use crate::boussinesq_solver::pre_built_components::
non_insulated_fluid_components::NonInsulatedFluidComponent;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::
heat_transfer_interaction_enums::HeatTransferInteractionType;


/// fluid mechanics bit 
/// calculate the fluid mechanics for the two branches in parallel
///
/// In actual fact though, it is just one branch and we are getting 
/// the mass flowrate through that branch,
/// can be used for
/// DHX + Heater branch (both branches form one loop)
pub fn get_mass_flowrate_across_two_branches(dhx_and_heater_branches: &FluidComponentSuperCollection) -> 
MassRate {
    let pressure_change_across_each_branch = 
        dhx_and_heater_branches.get_pressure_change(MassRate::ZERO);

    let mass_flowrate_across_each_branch: Vec<MassRate> = 
        dhx_and_heater_branches.
        get_mass_flowrate_across_each_parallel_branch(
            pressure_change_across_each_branch
        );

    let mut mass_flowrate: MassRate = 
        mass_flowrate_across_each_branch[0];


    // get absolute value
    mass_flowrate = mass_flowrate.abs();

    mass_flowrate

}

/// fluid mechanics calcs, specific to the primary (DHX plus Heater branch) loop
/// note that this only works if the components are correct
/// obtains mass flowrate across the primary (DHX plus Heater branch) loop
/// for hydrostatic pressure, all component angles are taken with 
/// reference to the branching point (around 5a)
///
/// note: only Flowmeter is  non insulated, all other components 
/// should be insulated
///
pub fn pri_loop_branches_fluid_mechanics_calc_mass_rate(
    pipe_4: &InsulatedFluidComponent,
    pipe_3: &InsulatedFluidComponent,
    pipe_2a: &InsulatedFluidComponent,
    static_mixer_10_label_2: &InsulatedFluidComponent,
    heater_top_head_1a: &InsulatedFluidComponent,
    heater_version1_1: &InsulatedFluidComponent,
    heater_bottom_head_1b: &InsulatedFluidComponent,
    pipe_18: &InsulatedFluidComponent,
    pipe_5a: &NonInsulatedFluidComponent,
    pipe_26: &InsulatedFluidComponent,
    pipe_25a: &NonInsulatedFluidComponent,
    static_mixer_21_label_25: &InsulatedFluidComponent,
    dhx_shell_side_pipe_24: &InsulatedFluidComponent,
    static_mixer_20_label_23: &InsulatedFluidComponent,
    pipe_23a: &InsulatedFluidComponent,
    pipe_22: &InsulatedFluidComponent,
    flowmeter_20_21a: &NonInsulatedFluidComponent,
    pipe_21: &InsulatedFluidComponent,
    pipe_20: &InsulatedFluidComponent,
    pipe_19: &InsulatedFluidComponent,
    pipe_17b: &InsulatedFluidComponent,
)-> MassRate {

    let mut heater_branch = 
        FluidComponentCollection::new_series_component_collection();

    heater_branch.clone_and_add_component(pipe_4);
    heater_branch.clone_and_add_component(pipe_3);
    heater_branch.clone_and_add_component(pipe_2a);
    heater_branch.clone_and_add_component(static_mixer_10_label_2);
    heater_branch.clone_and_add_component(heater_top_head_1a);
    heater_branch.clone_and_add_component(heater_version1_1);
    heater_branch.clone_and_add_component(heater_bottom_head_1b);
    heater_branch.clone_and_add_component(pipe_18);


    let mut dhx_branch = 
        FluidComponentCollection::new_series_component_collection();

    dhx_branch.clone_and_add_component(pipe_5a);
    dhx_branch.clone_and_add_component(pipe_26);
    dhx_branch.clone_and_add_component(pipe_25a);
    dhx_branch.clone_and_add_component(static_mixer_21_label_25);
    dhx_branch.clone_and_add_component(dhx_shell_side_pipe_24);
    dhx_branch.clone_and_add_component(static_mixer_20_label_23);
    dhx_branch.clone_and_add_component(pipe_23a);
    dhx_branch.clone_and_add_component(pipe_22);
    dhx_branch.clone_and_add_component(flowmeter_20_21a);
    dhx_branch.clone_and_add_component(pipe_21);
    dhx_branch.clone_and_add_component(pipe_20);
    dhx_branch.clone_and_add_component(pipe_19);
    dhx_branch.clone_and_add_component(pipe_17b);

    let mut pri_loop_branches = 
        FluidComponentSuperCollection::default();

    pri_loop_branches.set_orientation_to_parallel();
    pri_loop_branches.fluid_component_super_vector.push(heater_branch);
    pri_loop_branches.fluid_component_super_vector.push(dhx_branch);

    let mass_rate = get_mass_flowrate_across_two_branches(&pri_loop_branches);

    mass_rate

}

/// now the heat transfer for the DRACS loop 
/// for a single timestep, given mass flowrate in a counter clockwise 
/// fashion in the DRACS
///
/// you also must specify the heat transfer coefficient to ambient 
/// which is assumed to be the same throughout the loop
///
/// TODO: change this to actual DHX and Heater 
pub fn calculate_pri_loop_dhx_heater_thermal_hydraulics(
    mass_flowrate_counter_clockwise: MassRate,
    heat_rate_through_dhx: Power,
    tchx_heat_transfer_coeff: HeatTransfer,
    average_temperature_for_density_calcs: ThermodynamicTemperature,
    timestep: Time,
    ambient_htc: HeatTransfer,
    pipe_4: &InsulatedFluidComponent,
    pipe_3: &InsulatedFluidComponent,
    pipe_2a: &InsulatedFluidComponent,
    static_mixer_10_label_2: &InsulatedFluidComponent,
    heater_top_head_1a: &InsulatedFluidComponent,
    heater_version1_1: &InsulatedFluidComponent,
    heater_bottom_head_1b: &InsulatedFluidComponent,
    pipe_18: &InsulatedFluidComponent,
    pipe_5a: &NonInsulatedFluidComponent,
    pipe_26: &InsulatedFluidComponent,
    pipe_25a: &NonInsulatedFluidComponent,
    static_mixer_21_label_25: &InsulatedFluidComponent,
    dhx_shell_side_pipe_24: &InsulatedFluidComponent,
    static_mixer_20_label_23: &InsulatedFluidComponent,
    pipe_23a: &InsulatedFluidComponent,
    pipe_22: &InsulatedFluidComponent,
    flowmeter_20_21a: &NonInsulatedFluidComponent,
    pipe_21: &InsulatedFluidComponent,
    pipe_20: &InsulatedFluidComponent,
    pipe_19: &InsulatedFluidComponent,
    pipe_17b: &InsulatedFluidComponent,
    ){

        // for an ideal situation, we have zero parasitic heat losses
        // therefore, for each component, except tchx, heat transfer 
        // coeff is zero

        let adiabatic_heat_transfer_coeff = ambient_htc;

        // create the heat transfer interaction 
        let advection_heat_transfer_interaction: HeatTransferInteractionType;

        // I'm going to create the advection interaction
        //
        // and probably for the sake of density calcs, I'll take the 
        // average density using DHX outlet and 
        // TCHX outlet temperatures, average them for the whole loop 
        // doesn't make much diff tho based on Boussinesq approximation

        let average_therminol_density = 
            LiquidMaterial::TherminolVP1.density(
                average_temperature_for_density_calcs).unwrap();

        advection_heat_transfer_interaction = 
            HeatTransferInteractionType::
            new_advection_interaction(mass_flowrate_counter_clockwise, 
                average_therminol_density, 
                average_therminol_density);

        // now, let's link the fluid arrays using advection 
        // (no conduction here axially between arrays)
        {
        }
        // set the relevant heat transfer coefficients 
        // all zero except for tchx
        {
        }
        // add lateral heat losses and power through dhx
        {
        }

        // now we should be ready to advance timestep
        {
        }

        // we do it in serial, so it keeps things simple 
        // now we are done
        //
        todo!()

}
