
use std::f64::consts::PI;

use uom::si::f64::*;
use uom::si::thermodynamic_temperature::kelvin;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::density::try_get_rho;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::specific_enthalpy::try_get_h;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::specific_enthalpy::try_get_temperature_from_h;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::thermal_diffusivity::try_get_alpha_thermal_diffusivity;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::Material;
use crate::boussinesq_solver::control_volume_dimensions::InnerDiameterThermalConduction;
use crate::boussinesq_solver::control_volume_dimensions::OuterDiameterThermalConduction;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_geometry::CylindricalAndSphericalSolidFluidArrangement;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::*;
use crate::boussinesq_solver::single_control_vol::SingleCVNode;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::{DataAdvection, HeatTransferInteractionType};

use uom::num_traits::Zero;

use super::bc_types::BCType;
use super::cv_types::CVType;
use super::HeatTransferEntity;
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
#[inline]
pub fn link_heat_transfer_entity(entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType)-> Result<(), ThermalHydraulicsLibError>{

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
        (
            HeatTransferEntity::ControlVolume(cv_type_1),
            HeatTransferEntity::ControlVolume(cv_type_2)
        ) => {
            calculate_control_volume_serial(
                cv_type_1, 
                cv_type_2, 
                interaction)
        },
        (
            HeatTransferEntity::BoundaryConditions(bc_type),
            HeatTransferEntity::ControlVolume(cv_type)
        ) =>
        {
            attach_boundary_condition_to_control_volume_back_serial(
                bc_type, cv_type, interaction)
        },
        (
            HeatTransferEntity::ControlVolume(cv_type),
            HeatTransferEntity::BoundaryConditions(bc_type)
        ) =>
        {
            attach_boundary_condition_to_control_volume_front_serial(
                cv_type, bc_type, interaction)
        },
        (
            HeatTransferEntity::BoundaryConditions(bc_type_1),
            HeatTransferEntity::BoundaryConditions(bc_type_2)
        ) =>
        {
            return Err( 
                ThermalHydraulicsLibError::NotImplementedForBoundaryConditions
                ( "interactions between two boundary conditions \n not implemented".to_string()));
        },
    };


    return heat_transfer_entity_result;
   
}


// the job of this function is to take in a control volume 
// and then mutate it by calculating its interaction
#[inline]
pub (crate) fn calculate_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<(),ThermalHydraulicsLibError>{

    // let me first match my control volumes to their various types
    // at each matching arm, use those as inputs

    let cv_result = match (control_vol_1, control_vol_2) {
        (
            CVType::SingleCV(single_cv_1), CVType::SingleCV(single_cv_2)
        ) =>
        {
            calculate_between_two_singular_cv_nodes(
                single_cv_1, 
                single_cv_2, 
                interaction)
        },
        // basically like this 
        //
        // (array_cv_back -----array_cv --- array_cv_front) --- (single_cv) 
        (CVType::SolidArrayCV(solid_array_cv), CVType::SingleCV(single_cv)) =>
        {
            solid_array_cv.link_single_cv_to_higher_side(
                single_cv,
                interaction)
        },
        // basically like this 
        //
        // (array_cv_back -----array_cv --- array_cv_front) --- (single_cv) 
        (CVType::FluidArrayCV(fluid_array_cv), CVType::SingleCV(single_cv)) =>
        {
            fluid_array_cv.link_single_cv_to_higher_side(
                single_cv,
                interaction)
        },
        // basically like this 
        //
        // (single_cv) --- (array_cv_back -----array_cv --- array_cv_front)
        (CVType::SingleCV(single_cv), CVType::SolidArrayCV(solid_array_cv)) =>
        {
            solid_array_cv.link_single_cv_to_lower_side(
                single_cv,
                interaction)
        },
        // basically like this 
        //
        // (single_cv) --- (array_cv_back -----array_cv --- array_cv_front)
        (CVType::SingleCV(single_cv), CVType::FluidArrayCV(fluid_array_cv)) =>
        {
            fluid_array_cv.link_single_cv_to_lower_side(
                single_cv,
                interaction)
        },
        // basically like this 
        //
        // (back --- cv_1 --- front) ---- (back --- cv_2 --- front)
        (CVType::SolidArrayCV(solid_array_cv_1), CVType::SolidArrayCV(solid_array_cv_2)) =>
        {
            solid_array_cv_1.link_solid_column_to_the_front_of_this_solid_column(
                solid_array_cv_2,
                interaction)
        },
        // basically like this 
        //
        // (back --- cv_1 --- front) ---- (back --- cv_2 --- front)
        (CVType::SolidArrayCV(solid_array_cv_1), CVType::FluidArrayCV(fluid_array_cv_2)) =>
        {
            solid_array_cv_1.link_fluid_array_to_the_front_of_this_solid_column(
                fluid_array_cv_2,
                interaction)
        },
        // basically like this 
        //
        // (back --- cv_1 --- front) ---- (back --- cv_2 --- front)
        (CVType::FluidArrayCV(fluid_array_cv_1), CVType::SolidArrayCV(solid_array_cv_2)) =>
        {
            fluid_array_cv_1.link_solid_column_to_the_front_of_this_fluid_array(
                solid_array_cv_2,
                interaction)
        },
        // basically like this 
        //
        // (back --- cv_1 --- front) ---- (back --- cv_2 --- front)
        (CVType::FluidArrayCV(fluid_array_cv_1), CVType::FluidArrayCV(fluid_array_cv_2)) =>
        {
            fluid_array_cv_1.link_fluid_array_to_the_front_of_this_fluid_array(
                fluid_array_cv_2,
                interaction)
        },
    };


    return cv_result;

}

/// which calls other functions depending on whether the 
/// heat transfer interaction is conductance based on advection based
#[inline]
pub fn calculate_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType)-> Result<(), ThermalHydraulicsLibError>{


    match interaction {
        HeatTransferInteractionType::UserSpecifiedThermalConductance(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::SingleCartesianThermalConductanceOneDimension(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCartesianThermalConductanceThreeDimension(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCartesianThermalConductance(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::DualCylindricalThermalConductance(_, _, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::CylindricalConductionConvectionLiquidOutside(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::CylindricalConductionConvectionLiquidInside(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatAddition => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCustomArea(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCylindricalOuterArea(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedHeatFluxCylindricalInnerArea(_, _) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },
        HeatTransferInteractionType::UserSpecifiedConvectionResistance(_) => {
            calculate_conductance_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                interaction)
        },

        HeatTransferInteractionType::Advection(advection_data) => {
            calculate_advection_interaction_between_two_singular_cv_nodes(
                single_cv_1,
                single_cv_2,
                advection_data)
        },
    }

}
#[inline]
pub fn calculate_advection_interaction_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    advection_data: DataAdvection)-> Result<(), ThermalHydraulicsLibError>{

    let mass_flow_from_cv_1_to_cv_2 = advection_data.mass_flowrate;

    // for this, quite straightforward, 
    // get both specific enthalpy of both cvs 
    // calculate power flow 
    // and then update the power vector 
    //

    let specific_enthalpy_cv1: AvailableEnergy = 
    single_cv_1.current_timestep_control_volume_specific_enthalpy;

    let specific_enthalpy_cv2: AvailableEnergy = 
    single_cv_2.current_timestep_control_volume_specific_enthalpy;

    // calculate heat rate 

    let heat_flowrate_from_cv_1_to_cv_2: Power 
    = advection_heat_rate(mass_flow_from_cv_1_to_cv_2,
        specific_enthalpy_cv1,
        specific_enthalpy_cv2,)?;

    // by default, cv 1 is on the left, cv2 is on the right 
    //

    single_cv_1.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_cv_1_to_cv_2);
    single_cv_2.rate_enthalpy_change_vector.
        push(heat_flowrate_from_cv_1_to_cv_2);

    // relevant timescale here is courant number
    //
    // the timescale can only be calculated after the mass flows 
    // in and out of the cv are sufficiently calculated
    // the only thing we can do here is push the mass flowrate 
    // into the individual mass flowrate vectors 
    //
    // by convention, mass flowrate goes out of cv1 and into cv2 
    // so positive mass flowrate here is positive for cv2 
    // and negative for cv1
    //
    // I'll need a density for the flow first

    let density_cv1 = advection_data.fluid_density_heat_transfer_entity_1;
    let density_cv2 = advection_data.fluid_density_heat_transfer_entity_2;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_cv_1_to_cv_2 > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from cv1 
        // to cv2 
        // then the density we use is cv1 

        volumetric_flowrate = mass_flow_from_cv_1_to_cv_2/density_cv1;

    } else {
        // if mass flowrate is positive, flow is moving from cv2
        // to cv1
        // then the density we use is cv2

        volumetric_flowrate = mass_flow_from_cv_1_to_cv_2/density_cv2;
    }

    // now that I've done the volume flowrate calculation, push the 
    // volumetric flowrate to each vector
    single_cv_1.volumetric_flowrate_vector.push(
        -volumetric_flowrate);
    single_cv_2.volumetric_flowrate_vector.push(
        volumetric_flowrate);

    

    // done! 
    Ok(())
}

#[inline]
pub fn calculate_conductance_interaction_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType)-> Result<(), ThermalHydraulicsLibError>{

    // let's get the two temperatures of the control volumes first
    // so let me get the enthalpies, and then their respective 
    // temperatures 

    let single_cv_1_enthalpy = single_cv_1.
        current_timestep_control_volume_specific_enthalpy;
    let single_cv_2_enthalpy = single_cv_2.
        current_timestep_control_volume_specific_enthalpy;

    // to get the temperatures, we'll need the material as well 
    let single_cv_1_material = single_cv_1.material_control_volume;
    let single_cv_2_material = single_cv_2.material_control_volume;

    // we'll also need to get their pressures 
    let single_cv_1_pressure = single_cv_1.pressure_control_volume;
    let single_cv_2_pressure = single_cv_2.pressure_control_volume;

    // we will now get their respective temperatures 
    //
    // (note, this is extremely computationally expensive as it 
    // is iterative in nature)
    //
    // two solutions here, 
    // one: store cv temperature in single cv, 
    // so it can be readily accessed
    //
    // two: cheaper method of getting t from h.
    

    let single_cv_1_temperature: ThermodynamicTemperature;
    let single_cv_2_temperature: ThermodynamicTemperature;

    let experimental_code = true; 
    if experimental_code {

        single_cv_1_temperature = single_cv_1.temperature;
        single_cv_2_temperature = single_cv_2.temperature;

    } else {

        // original code
        single_cv_1_temperature = try_get_temperature_from_h(
            single_cv_1_material, 
            single_cv_1_enthalpy, 
            single_cv_1_pressure)?;
        single_cv_2_temperature = try_get_temperature_from_h(
            single_cv_2_material, 
            single_cv_2_enthalpy, 
            single_cv_2_pressure)?;

    }

    // now that we got their respective temperatures we can calculate 
    // the thermal conductance between them
    //
    // for conduction for instance, q = kA dT/dx 
    // conductance is watts per kelvin or 
    // q = (kA)/dx * dT
    // conductance here is kA/dx
    // thermal resistance is 1/conductance
    //
    // for convection, we get: 
    // q = hA (Delta T)
    // hA becomes the thermal conductance
    //
    // If we denote thermal conductance as Htc
    // 
    // Then a general formula for heat flowing from 
    // temperature T_1 to T_2 is 
    //
    // T_1 --> q --> T_2 
    //
    // q = - Htc (T_2 - T_1)

    // 
    let thermal_conductance = get_thermal_conductance_based_on_interaction(
        single_cv_1_temperature, 
        single_cv_2_temperature,
        single_cv_1_pressure, 
        single_cv_2_pressure, 
        interaction)?;

    // suppose now we have thermal conductance, we can now obtain the 
    // power flow
    //

    let cv_2_temp_minus_cv_1_temp_kelvin: f64 = 
        single_cv_2_temperature.get::<kelvin>() - 
        single_cv_1_temperature.get::<kelvin>();

    let cv_2_temp_minus_cv_1: TemperatureInterval = 
    TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(
        cv_2_temp_minus_cv_1_temp_kelvin);

    let heat_flowrate_from_cv_1_to_cv_2: Power = 
    - thermal_conductance * cv_2_temp_minus_cv_1;

    // now, we add a heat loss term to cv_1 
    // and a heat gain term to cv_2 
    //
    // using timestep
    // the signs should cancel out

    single_cv_1.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_cv_1_to_cv_2);
    single_cv_2.rate_enthalpy_change_vector.
        push(heat_flowrate_from_cv_1_to_cv_2);


    // for solids mesh fourier number need only 
    // be done once, not every time 
    // an interaction is formed 
    //
    // probably the cell stability fourier number will be done in the 
    // constructor. however, with convection, the time scale must be 
    // recalculated at every time step. so it really depends whether 
    // it's solid or fluid control volume
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    //
    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    //
    // but the liquid timescales are calculated at the cv level only 
    // after all volumetric flowrates are calculated
    //
    // so don't really need any new timescale calculations

    return Ok(());

}
/// this function calculates relevant timescales when linking 
/// two heat transfer entities
pub fn calculate_timescales_for_heat_transfer_entity(
    entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType)-> Result<Time, ThermalHydraulicsLibError>{

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
            return Err( 
                ThermalHydraulicsLibError::NotImplementedForBoundaryConditions
                ( "interactions between two boundary conditions \n not implemented".to_string()));
        }
    };


    return heat_transfer_entity_timestep_result;
   
}

/// calculates timestep for control voluem and boundary condition 
/// in serial
/// for arrayCVs, the boundary condition is attached to the front 
/// of the control volume
///
/// (back --- cv_1 --- front) ---- (boundary condition)
pub (crate) fn calculate_timestep_control_volume_front_to_boundary_condition_serial(
    control_vol: &mut CVType,
    boundary_condition: &mut BCType,
    interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError> {


    let cv_bc_result = match (control_vol, boundary_condition) {
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(_heat_rate)) 
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(_heat_flux))
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedTemperature(_bc_temperature))
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (&mut solid_array_cv.front_single_cv, interaction)
        }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut solid_array_cv.front_single_cv, interaction)
        }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut solid_array_cv.front_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (&mut fluid_array_cv.front_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut fluid_array_cv.front_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut fluid_array_cv.front_single_cv, interaction)
        }
        ,
    };

    return cv_bc_result;
}

/// the job of this function is to calculate time step 
/// when two control volumes interact, the control volumes can be 
/// single control volumes or array control volumes
#[inline]
pub (crate) fn calculate_timestep_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError>{

    // let me first match my control volumes to their various types
    // at each matching arm, use those as inputs

    let cv_result = match (control_vol_1, control_vol_2) {
        (CVType::SingleCV(single_cv_1), CVType::SingleCV(single_cv_2)) =>
        {
            single_cv_1.calculate_mesh_stability_timestep_for_two_single_cv_nodes(
                single_cv_2, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv), CVType::SingleCV(single_cv_node)) =>
        {
            fluid_array_cv.calculate_timestep_for_single_cv_to_front_of_array_cv(
                    single_cv_node,
                    interaction)
        }
        ,
        (CVType::SolidArrayCV(solid_array_cv), CVType::SingleCV(single_cv_node)) =>
        {
            solid_array_cv.calculate_timestep_for_single_cv_to_front_of_array_cv(
                    single_cv_node,
                    interaction)
        }
        ,
        (CVType::SingleCV(single_cv_node), CVType::SolidArrayCV(solid_array_cv)) =>
        {
            solid_array_cv.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    single_cv_node,
                    interaction)
        },
        (CVType::SingleCV(single_cv_node), CVType::FluidArrayCV(fluid_array_cv)) =>
        {
            fluid_array_cv.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    single_cv_node,
                    interaction)
        },
        (CVType::FluidArrayCV(fluid_array_cv_1), CVType::FluidArrayCV(fluid_array_cv_2)) =>
        {
            fluid_array_cv_2.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    &mut fluid_array_cv_1.front_single_cv,
                    interaction)
        },
        (CVType::SolidArrayCV(solid_array_cv_1), CVType::SolidArrayCV(solid_array_cv_2)) =>
        {
            solid_array_cv_2.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    &mut solid_array_cv_1.front_single_cv,
                    interaction)
        },
        (CVType::FluidArrayCV(fluid_array_cv_1), CVType::SolidArrayCV(solid_array_cv_2)) =>
        {
            solid_array_cv_2.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    &mut fluid_array_cv_1.front_single_cv,
                    interaction)
        },
        (CVType::SolidArrayCV(solid_array_cv_1), CVType::FluidArrayCV(fluid_array_cv_2)) =>
        {
            fluid_array_cv_2.
                calculate_timestep_for_single_cv_to_back_of_array_cv(
                    &mut solid_array_cv_1.front_single_cv,
                    interaction)
        },
    };


    return cv_result;

}

/// calculates timestep for control voluem and boundary condition 
/// in serial
/// for arrayCVs, the boundary condition is attached to the back
/// of the control volume
///
/// (boundary condition) ----- (back --- cv_1 --- front) 
///
pub (crate) fn calculate_timestep_control_volume_back_to_boundary_condition_serial(
    boundary_condition: &mut BCType,
    control_vol: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError> {


    let cv_bc_result = match (control_vol, boundary_condition) {
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(_heat_rate)) 
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(_heat_flux))
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedTemperature(_bc_temperature))
            => {
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (single_cv, interaction)
            }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (&mut solid_array_cv.back_single_cv, interaction)
        }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut solid_array_cv.back_single_cv, interaction)
        }
        ,
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut solid_array_cv.back_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                    (&mut fluid_array_cv.back_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut fluid_array_cv.back_single_cv, interaction)
        }
        ,
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            // match the cv to the cv type, extract the front boundary 
            // condition 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc
                (&mut fluid_array_cv.back_single_cv, interaction)
        }
        ,
    };

    return cv_bc_result;
}

/// this is thermal conductance function based on interaction type 
/// may want to move the calculation bits to the calculation module in 
/// future
/// it calculates thermal conductance based on the supplied enum
///
/// TODO: probably want to test this function out
/// 
pub (in crate) 
fn get_thermal_conductance_based_on_interaction(
    temperature_1: ThermodynamicTemperature,
    temperature_2: ThermodynamicTemperature,
    pressure_1: Pressure,
    pressure_2: Pressure,
    interaction: HeatTransferInteractionType) 
-> Result<ThermalConductance, ThermalHydraulicsLibError> 
{

    let conductance: ThermalConductance = match 
        interaction {
            HeatTransferInteractionType::UserSpecifiedThermalConductance(
                user_specified_conductance) => user_specified_conductance,
            HeatTransferInteractionType
                ::SingleCartesianThermalConductanceOneDimension(
                material,thickness) => get_conductance_single_cartesian_one_dimension(
                    material,
                    temperature_1, 
                    temperature_2, 
                    pressure_1, 
                    pressure_2, 
                    thickness)?,
            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(
                (solid_material, shell_thickness,
                solid_temperature, solid_pressure), 
                (h, inner_diameter, cylinder_length)) => {

                    let id: Length = inner_diameter.clone().into();
                    let thicnkess: Length = shell_thickness.clone().into();

                    let od: Length = id+thicnkess;

                    let outer_diameter: OuterDiameterThermalConduction = 
                    OuterDiameterThermalConduction::from(od);

                    // after all the typing conversion, we can 
                    // get our conductance
                    get_conductance_single_cylindrical_radial_solid_liquid(
                        solid_material,
                        solid_temperature,
                        solid_pressure,
                        h,
                        inner_diameter,
                        outer_diameter,
                        cylinder_length,
                        CylindricalAndSphericalSolidFluidArrangement::
                        FluidOnInnerSurfaceOfSolidShell ,
                    )?
                },

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidOutside(
                (solid_material, shell_thickness,
                solid_temperature, solid_pressure), 
                (h, outer_diameter, cylinder_length)) => {

                    let od: Length = outer_diameter.clone().into();
                    let thicnkess: Length = shell_thickness.clone().into();

                    let id: Length = od - thicnkess;

                    let inner_diameter: InnerDiameterThermalConduction = 
                    InnerDiameterThermalConduction::from(id);

                    // after all the typing conversion, we can 
                    // get our conductance
                    get_conductance_single_cylindrical_radial_solid_liquid(
                        solid_material,
                        solid_temperature,
                        solid_pressure,
                        h,
                        inner_diameter,
                        outer_diameter,
                        cylinder_length,
                        CylindricalAndSphericalSolidFluidArrangement::
                        FluidOnInnerSurfaceOfSolidShell ,
                    )?
                },
            // note: actually function signatures are a little more 
            // friendly to use than packing enums with lots of stuff 
            // so may change stuffing enums with tuples to stuffing 
            // enums with a single struct
            HeatTransferInteractionType
                ::DualCylindricalThermalConductance(
                (inner_material,inner_shell_thickness),
                (outer_material,outer_shell_thickness),
                (inner_diameter,
                outer_diameter,
                cylinder_length)
            ) => {
                    // first, want to check if inner_diameter + 
                    // shell thicknesses is outer diameter 

                    let expected_outer_diameter: Length;
                    let id: Length = inner_diameter.into();
                    let inner_thickness: Length =  inner_shell_thickness.into();
                    let outer_thickness: Length =  outer_shell_thickness.into();

                    expected_outer_diameter = 
                        id + inner_thickness + outer_thickness;

                    let od: Length = outer_diameter.into();

                    // inner diameter and outer diameter values must be 
                    // equal to within 1 nanometer 1e-9 m
                    if (od.value - expected_outer_diameter.value).abs() > 1e-9
                    {

                        let mut error_str: String = "the inner diameter 
                            plus shell thicknesses do not equate 
                            to outer diameter".to_string();

                        error_str += "supplied outer diameter (m):";
                        error_str += &od.value.to_string();
                        error_str += "expected outer diameter (m):";
                        error_str += &expected_outer_diameter.value.to_string();


                        return Err(ThermalHydraulicsLibError::GenericStringError(error_str));
                    }

                    get_conductance_cylindrical_radial_two_materials(
                        inner_material,
                        outer_material,
                        temperature_1, //convention, 1 is inner shell
                        temperature_2, // convention 2, is outer shell
                        pressure_1,
                        pressure_2,
                        inner_diameter,
                        inner_shell_thickness,
                        outer_shell_thickness,
                        cylinder_length,
                    )?
                },
            HeatTransferInteractionType::UserSpecifiedHeatAddition  
                => {
                    println!("interaction type needs to be \n 
                        thermal conductance");
                    return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
                },

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCustomArea(_) => {

                    println!("interaction type needs to be \n 
                        thermal conductance");
                    return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(_,_) => {

                    println!("interaction type needs to be \n 
                        thermal conductance");
                    return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalInnerArea(_,_) => {

                    println!("interaction type needs to be \n 
                        thermal conductance");
                    return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);

                },
            HeatTransferInteractionType::
                DualCartesianThermalConductance(
                (material_1, thickness_1),
                (material_2,thickness_2)) => { 
                    
                    let conductnace_layer_1: ThermalConductance 
                    = get_conductance_single_cartesian_one_dimension(
                        material_1,
                        temperature_1, 
                        temperature_1, 
                        pressure_1, 
                        pressure_1, 
                        thickness_1)?;

                    let conductnace_layer_2: ThermalConductance 
                    = get_conductance_single_cartesian_one_dimension(
                        material_2,
                        temperature_2, 
                        temperature_2, 
                        pressure_2, 
                        pressure_2, 
                        thickness_2)?;

                    let overall_resistance = 
                    1.0/conductnace_layer_2 
                    + 1.0/conductnace_layer_1;

                    // return the conductance or resistnace inverse

                    1.0/overall_resistance
            },
            HeatTransferInteractionType::
                DualCartesianThermalConductanceThreeDimension(
                data_dual_cartesian_conduction) 
                => {

                    let material_1 = 
                    data_dual_cartesian_conduction .material_1;

                    let material_2 = 
                    data_dual_cartesian_conduction .material_2;

                    let thickness_1 = 
                    data_dual_cartesian_conduction .thickness_1;

                    let thickness_2 = 
                    data_dual_cartesian_conduction .thickness_2;

                    let xs_area = 
                    data_dual_cartesian_conduction .xs_area;

                    get_conductance_dual_cartesian_three_dimensions(
                        material_1, 
                        material_2, 
                        temperature_1, 
                        temperature_2, 
                        pressure_1, 
                        pressure_2, 
                        xs_area, 
                        thickness_1,
                        thickness_2)?
                },

            HeatTransferInteractionType::
                UserSpecifiedConvectionResistance(
                data_convection_resistance) 
                => {

                    let heat_transfer_coeff: HeatTransfer = 
                    data_convection_resistance.heat_transfer_coeff;
                    let surf_area: Area = 
                    data_convection_resistance.surf_area.into();

                    heat_transfer_coeff * surf_area
                },

            HeatTransferInteractionType::Advection(_) => {
                println!("advection interaction types \n 
                do not correspond to conductance");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            },

        };

    return Ok(conductance);
}

/// the job of this function is to handle interactions between 
/// a control_volume and a boundary condition
///
/// here we have to handle some BC types, 
/// these are the most basic:
/// 1. constant heat flux
/// 2. constant temperature
/// 3. constant heat addition
///
/// For each case, there should be a function
/// to handle each case
/// 
/// for arrayCVs, the boundary condition is attached to the front 
/// of the control volume
///
/// (back --- cv_1 --- front) ---- (boundary condition)
#[inline]
pub (crate) fn attach_boundary_condition_to_control_volume_front_serial(
    control_vol: &mut CVType,
    boundary_condition: &mut BCType,
    interaction: HeatTransferInteractionType) -> Result<(),ThermalHydraulicsLibError> {

    let cv_bc_result = match (control_vol, boundary_condition) {
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(heat_rate)) 
            => calculate_constant_heat_addition_front_single_cv_back(
                single_cv, *heat_rate, interaction),
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(heat_flux))
            => calculate_constant_heat_flux_front_single_cv_back(
                single_cv, *heat_flux, interaction),
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedTemperature(bc_temperature))
            => calculate_constant_temperature_front_single_cv_back(
                single_cv, *bc_temperature, interaction),
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            solid_array_cv.link_heat_flux_bc_to_front_of_this_cv(
                *heat_flux,
                interaction)
        },
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            solid_array_cv.link_heat_addition_to_front_of_this_cv(
                *heat_rate,
                interaction)
        },
        (CVType::SolidArrayCV(solid_array_cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            solid_array_cv.link_constant_temperature_to_front_of_this_cv(
                *bc_temperature,
                interaction)
        },
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            fluid_array_cv.link_heat_flux_bc_to_front_of_this_cv(
                *heat_flux,
                interaction)
        },
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            fluid_array_cv.link_heat_addition_to_front_of_this_cv(
                *heat_rate,
                interaction)
        },
        (CVType::FluidArrayCV(fluid_array_cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            fluid_array_cv.link_constant_temperature_to_front_of_this_cv(
                *bc_temperature,
                interaction)
        },
    };

    return cv_bc_result;
}

/// for connecting a bc to cv where 
///
/// (cv) ---------- (constant temperature bc)
#[inline]
pub fn calculate_constant_temperature_front_single_cv_back(
    control_vol: &mut SingleCVNode,
    boundary_condition_temperature: ThermodynamicTemperature,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError> {

    match interaction {
        HeatTransferInteractionType::Advection(
            advection_dataset) => {

            // I'm mapping my own error to string, so off
            calculate_bc_front_cv_back_advection(
                control_vol,
                advection_dataset)?;
            return Ok(());
        },
        _ => (),
    }
    // if anything else, use conductance

    control_vol.calculate_single_cv_node_constant_temperature_conductance(
        boundary_condition_temperature,
        interaction)?;
    return Ok(());
}


/// calculates the interaction between a heat flux BC and 
/// a control volume 
///
/// (single cv) ------------------ (heat flux bc)
///
/// the heat addition is at the front, the cv is at the back
pub fn calculate_constant_heat_flux_front_single_cv_back(
    control_vol: &mut SingleCVNode,
    heat_flux_into_control_vol: HeatFluxDensity,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError> {

    // first, obtain a heat transfer area from the constant heat flux 
    // BC
    let heat_transfer_area: Area = match interaction{
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCartesianThermalConductance(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(_, _, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,
        // these interaction types are acceptable
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => area,

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(
                    cylinder_length, od) => {
                    let od: Length = od.into();
                    let cylinder_length: Length  = cylinder_length.into();

                    let area = PI * od * cylinder_length;
                    area
                }
        ,

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(
                cylinder_length, id) => {
                let id: Length = id.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * id * cylinder_length;
                area

            }
        ,

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,
        HeatTransferInteractionType::
            Advection(advection_data) => {

                calculate_bc_front_cv_back_advection(
                    control_vol,
                    advection_data)?;

                return Ok(());
            }
        ,
    };

    let heat_flowrate_into_control_vol: Power = 
        heat_flux_into_control_vol * heat_transfer_area;

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_into_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change

    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    //
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Material::Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;
            ()
        },
        Material::Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }
    
    return Ok(());


}

/// calculates the interaction between a heat flux BC and 
/// a control volume 
///
/// (heat flux bc) ------------------ (single cv)
///
/// the cv is at the front 
/// heat addition is at the back
pub fn calculate_single_cv_front_heat_flux_back(
    heat_flux_into_control_vol: HeatFluxDensity,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError> {

    // first, obtain a heat transfer area from the constant heat flux 
    // BC
    let heat_transfer_area: Area = match interaction{
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCartesianThermalConductance(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(_, _, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(_, _) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,

        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        ,
        // these interaction types are acceptable
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => area,

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(
            cylinder_length, od) => {
                let od: Length = od.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * od * cylinder_length;
                area
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(
            cylinder_length, id) => {
                let id: Length = id.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * id * cylinder_length;
                area

            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => 
            {
                println!("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar");
                return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
            }
        HeatTransferInteractionType::
            Advection(advection_data) => 
            {
                calculate_cv_front_bc_back_advection(
                    control_vol,
                    advection_data)?;

                return Ok(());
            },
    };

    let heat_flowrate_into_control_vol: Power = 
    heat_flux_into_control_vol * heat_transfer_area;

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_into_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change

    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    //
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Material::Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;
            ()
        },
        Material::Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }
    
    return Ok(());


}

/// calculates the interaction between a heat addition BC and 
/// a control volume 
///
/// (single cv) ------------------ (heat addition bc)
///
/// the heat addition is at the front, the cv is at the back
#[inline]
pub fn calculate_constant_heat_addition_front_single_cv_back(
    control_vol: &mut SingleCVNode,
    heat_added_to_control_vol: Power,
    interaction: HeatTransferInteractionType
    ) -> Result<(), ThermalHydraulicsLibError> {

    // ensure that the interaction is UserSpecifiedHeatAddition
    // or advection
    // otherwise, return error 

    match interaction {

        HeatTransferInteractionType::UserSpecifiedHeatAddition => {
            // return a void value, that would be dropped 
            // instantly
            //
            // it pretty much has the same meaning as break
            ()
        },
        HeatTransferInteractionType::Advection(advection_data) => {
            calculate_bc_front_cv_back_advection(
                control_vol,
                advection_data)?;

            return Ok(());
            
        },


        _ => {
            println!("you need to specify that the interaction type \n 
            is UserSpecifiedHeatAddition");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);

        },
    };


    control_vol.rate_enthalpy_change_vector.
        push(heat_added_to_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change
    //
    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    


    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Material::Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;

            ()
        },
        Material::Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }

    return Ok(());
}

/// calculates the interaction between a heat addition BC and 
/// a control volume 
///
/// (heat addition) ------------------ (single cv)
///
/// the heat addition is at the front, the cv is at the back
#[inline]
pub fn calculate_single_cv_front_constant_heat_addition_back(
    heat_added_to_control_vol: Power,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType
    ) -> Result<(), ThermalHydraulicsLibError> {

    // ensure that the interaction is UserSpecifiedHeatAddition
    // or advection
    // otherwise, return error 

    match interaction {

        HeatTransferInteractionType::UserSpecifiedHeatAddition => {
            // return a void value, that would be dropped 
            // instantly
            //
            // it pretty much has the same meaning as break
            ()
        },
        HeatTransferInteractionType::Advection(advection_data) => {
            calculate_cv_front_bc_back_advection(
                control_vol,
                advection_data)?;

            return Ok(());

        },

        _ => {
            println!("you need to specify that the interaction type \n 
            is UserSpecifiedHeatAddition");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);

        },
    };


    control_vol.rate_enthalpy_change_vector.
        push(heat_added_to_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change
    //
    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    


    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Material::Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;

            ()
        },
        Material::Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }

    return Ok(());
}

pub fn calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError> {

    // here we have timestep based on the generic lengthscale of the 
    // control volume 
    let mut cv_timestep:Time = 
    control_vol.calculate_conduction_timestep()?;

    // we may have other time scales based on differing length scales 
    // of the control volume 
    //
    // so we will need to calculate time scales based on these other 
    // length scales and then calculate each of their own time scales.
    // if shorter, then we need to append it to the respective control 
    // volumes

    let cv_material = control_vol.material_control_volume.clone();
    let cv_pressure = control_vol.pressure_control_volume.clone();
    let cv_temperature = control_vol.temperature;

    
    let cv_alpha: DiffusionCoefficient = 
    try_get_alpha_thermal_diffusivity(cv_material,
        cv_temperature,
        cv_pressure)?;


    let max_mesh_fourier_number: f64 = 0.25;


    match interaction {
        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => {

                // do nothing

                ()
            },
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => {

                // if a conductance is specified, don't 
                // do anything

            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => {
                // when a normal area is given,
                // we can calculate volume to area ratio 

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }

            },

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(
            material, x_thickness) => {

                // the given material here overrides the normal 
                // material 
                let cv_alpha: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material,
                    cv_temperature,
                    cv_pressure)?;

                // if we connect the cv to a boundary condition,
                // then the length provided here is what we need 
                // to bother with 

                let lengthscale: Length = x_thickness.into();

                let time_step_max_based_on_x_thickness: Time 
                = max_mesh_fourier_number *
                lengthscale * 
                lengthscale / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_x_thickness {
                    cv_timestep = time_step_max_based_on_x_thickness;
                }
            },

        HeatTransferInteractionType::
            DualCartesianThermalConductance(
            (material_1, length_1), 
            (material_2,length_2)) => {
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;

                let length_1: Length = length_1.into();
                let length_2: Length = length_2.into();

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                length_2 *
                length_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()
            },

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(
            (material_1, radius_1), 
            (material_2,radius_2), 
            _) => {
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more radiuss
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;

                let radius_1: Length = radius_1.into();
                let radius_2: Length = radius_2.into();

                let timestep_1: Time = max_mesh_fourier_number * 
                radius_1 *
                radius_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                radius_2 *
                radius_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()
            },


        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(
            data_dual_cartesian_conduction_data) => {

                let material_1 = data_dual_cartesian_conduction_data.
                    material_1.clone();

                let material_2 = data_dual_cartesian_conduction_data.
                    material_2.clone();


                let length_1 : Length = data_dual_cartesian_conduction_data.
                    thickness_1.clone().into();

                let length_2 : Length = data_dual_cartesian_conduction_data.
                    thickness_2.clone().into();
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more radiuss
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;


                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                length_2 *
                length_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()

            },

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
            (material,radius,
            temperature,pressure),_) => {

                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cylindrical thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep

                let alpha_1: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material,
                    temperature,
                    pressure)?;

                let length_1: Length =  radius.into();

                let length_1 = length_1*0.5;

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }
                ()
                // if the control volume is fluid, we will need 
                // to introduce another time scale

            },

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(
            (material,radius,
            temperature,pressure),_) => {

                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cylindrical thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep

                let alpha_1: DiffusionCoefficient = 
                try_get_alpha_thermal_diffusivity(material,
                    temperature,
                    pressure)?;

                let length_1: Length =  radius.into();

                let length_1 = length_1*0.5;

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }
                ()
                // if the control volume is fluid, we will need 
                // to introduce another time scale

            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(l, od) => {

                // this is treated like a custom area kind of thing 
                // so we calculate the area first

                let cylinder_length: Length = l.into();
                let outer_diameter: Length = od.into();

                let area: Area = PI * outer_diameter * cylinder_length;

                // and then do the boilerplate code

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(l, id) => {

                // this is treated like a custom area kind of thing 
                // so we calculate the area first

                let cylinder_length: Length = l.into();
                let inner_diameter: Length = id.into();

                let area: Area = PI * inner_diameter * cylinder_length;

                // and then do the boilerplate code

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }

                
            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => {

                // if a resistance is specified, don't 
                // do anything

                ()
            },
        HeatTransferInteractionType::Advection(_) => {
            // advection has nothing to do with conduction timestep 
            // do nothing

            ()
        },

    }


    control_vol.max_timestep_vector.push(cv_timestep);

    return Ok(cv_timestep);

}

/// for advection calculations with heat flux or heat addition BC,
/// the temperature of flows flowing in and out of the BC will be 
/// determined by that of the control volume
///
/// it will be the same temperature as that of the control volume 
/// at that current timestep
///
/// this will be quite similar to how OpenFOAM treats inflows and outflows 
/// at zero gradient BCs
/// 
#[inline]
pub (crate) fn calculate_cv_front_bc_back_advection(
    control_vol: &mut SingleCVNode,
    advection_data: DataAdvection
) -> Result<(), ThermalHydraulicsLibError>{

    let mass_flow_from_bc_to_cv = advection_data.mass_flowrate;

    let specific_enthalpy_cv: AvailableEnergy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    // for the constant temperature bc, 
    // I'll just assume the fluid flowing from the bc is the same as the 
    // fluid in the cv
    //
    // and the pressure is same as the cv

    let control_vol_material = control_vol.material_control_volume;
    let control_vol_pressure = control_vol.pressure_control_volume;
    let control_vol_temperature = control_vol.get_temperature_from_enthalpy_and_set()?;

    let specific_enthalpy_bc_zero_gradient: AvailableEnergy = try_get_h(
        control_vol_material,
        control_vol_temperature,
        control_vol_pressure
    )?;
    // calculate heat rate 

    let heat_flowrate_from_bc_to_cv: Power 
    = advection_heat_rate(mass_flow_from_bc_to_cv,
        specific_enthalpy_bc_zero_gradient,
        specific_enthalpy_cv,)?;

    // push to cv
    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_from_bc_to_cv);


    let density_cv = advection_data.fluid_density_heat_transfer_entity_2;

    // we also need to calculate bc density 
    // now this doesn't quite work well for compressible flow but for 
    // now I'll let it be

    let density_bc: MassDensity = try_get_rho(
        control_vol_material,
        control_vol_temperature,
        control_vol_pressure
    )?;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_bc_to_cv > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from bc 
        // to cv 
        // then the density we use is bc 

        volumetric_flowrate = mass_flow_from_bc_to_cv/density_bc;

    } else {
        // if mass flowrate is positive, flow is moving from cv
        // to bc
        // then the density we use is cv

        volumetric_flowrate = mass_flow_from_bc_to_cv/density_cv;
    }


    // for courant number
    control_vol.volumetric_flowrate_vector.push(
        volumetric_flowrate);


    Ok(())
}
/// for advection calculations with heat flux or heat addition BC,
/// the temperature of flows flowing in and out of the BC will be 
/// determined by that of the control volume
///
/// it will be the same temperature as that of the control volume 
/// at that current timestep
/// 
/// this will be quite similar to how OpenFOAM treats inflows and outflows 
/// at zero gradient BCs
///
#[inline]
pub (crate) fn calculate_bc_front_cv_back_advection(
    control_vol: &mut SingleCVNode,
    advection_data: DataAdvection
) -> Result<(), ThermalHydraulicsLibError>{

    let mass_flow_from_cv_to_bc = advection_data.mass_flowrate;

    let specific_enthalpy_cv: AvailableEnergy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    // for the constant temperature bc, 
    // I'll just assume the fluid flowing from the bc is the same as the 
    // fluid in the cv
    //
    // and the pressure is same as the cv

    let control_vol_material = control_vol.material_control_volume;
    let control_vol_pressure = control_vol.pressure_control_volume;
    let control_vol_temperature = control_vol.get_temperature_from_enthalpy_and_set()?;

    let specific_enthalpy_bc_zero_gradient: AvailableEnergy = try_get_h(
        control_vol_material,
        control_vol_temperature,
        control_vol_pressure
    )?;
    // calculate heat rate 

    let heat_flowrate_from_bc_to_cv: Power 
    = advection_heat_rate(mass_flow_from_cv_to_bc,
        specific_enthalpy_cv,
        specific_enthalpy_bc_zero_gradient,)?;

    // push to cv
    control_vol.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_bc_to_cv);


    let density_cv = advection_data.fluid_density_heat_transfer_entity_2;

    // we also need to calculate bc density 
    // now this doesn't quite work well for compressible flow but for 
    // now I'll let it be

    let density_bc: MassDensity = try_get_rho(
        control_vol_material,
        control_vol_temperature,
        control_vol_pressure
    )?;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_cv_to_bc > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from bc 
        // to cv 
        // then the density we use is bc 

        volumetric_flowrate = mass_flow_from_cv_to_bc/density_bc;

    } else {
        // if mass flowrate is positive, flow is moving from cv
        // to bc
        // then the density we use is cv

        volumetric_flowrate = mass_flow_from_cv_to_bc/density_cv;
    }


    // for courant number
    control_vol.volumetric_flowrate_vector.push(
        -volumetric_flowrate);


    Ok(())
}

/// the job of this function is to handle interactions between 
/// a control_volume and a boundary condition
///
/// here we have to handle some BC types, 
/// these are the most basic:
/// 1. constant heat flux
/// 2. constant temperature
/// 3. constant heat addition
///
/// For each case, there should be a function
/// to handle each case,
/// for single CVs, it does pretty much the same job
/// 
/// for arrayCVs, the boundary condition is attached to the front 
/// of the control volume
///
/// (boundary condition) ---- (back --- cv_1 --- front) 
pub (crate) fn attach_boundary_condition_to_control_volume_back_serial(
    boundary_condition: &mut BCType,
    control_vol: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<(),ThermalHydraulicsLibError> {

    let cv_bc_result = match (control_vol, boundary_condition) {
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(heat_rate)) 
            => calculate_single_cv_front_constant_heat_addition_back(
                *heat_rate, single_cv, interaction),
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(heat_flux))
            => calculate_single_cv_front_heat_flux_back(
                *heat_flux,single_cv, interaction),
        (CVType::SingleCV(single_cv), BCType::UserSpecifiedTemperature(bc_temperature))
            => calculate_single_cv_node_front_constant_temperature_back(
                *bc_temperature, single_cv, interaction),
        (CVType::FluidArrayCV(cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            cv.link_heat_flux_bc_to_back_of_this_cv(
                *heat_flux,
                interaction)

        },
        (CVType::FluidArrayCV(cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            cv.link_heat_addition_to_back_of_this_cv(
                *heat_rate,
                interaction)
        },
        (CVType::FluidArrayCV(cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            cv.link_constant_temperature_to_back_of_this_cv(
                *bc_temperature,
                interaction)
        },
        (CVType::SolidArrayCV(cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            cv.link_heat_flux_bc_to_back_of_this_cv(
                *heat_flux,
                interaction)

        },
        (CVType::SolidArrayCV(cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            cv.link_heat_addition_to_back_of_this_cv(
                *heat_rate,
                interaction)
        },
        (CVType::SolidArrayCV(cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            cv.link_constant_temperature_to_back_of_this_cv(
                *bc_temperature,
                interaction)
        },
    };

    return cv_bc_result;
}



/// calculates interaction for a single cv and a constant temperature BC
#[inline]
pub fn calculate_single_cv_node_front_constant_temperature_back(
    boundary_condition_temperature: ThermodynamicTemperature,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError> {
    // this code is pretty crappy but I'll match advection first

    match interaction {
        HeatTransferInteractionType::Advection(
        advection_dataset) => {

                // I'm mapping my own error to string, so off
                calculate_cv_front_bc_back_advection(
                    control_vol,
                    advection_dataset)?;
                return Ok(());
            },
        _ => (),
    }

    // if anything else, use conductance

    control_vol.calculate_single_cv_node_constant_temperature_conductance(
        boundary_condition_temperature,
        interaction)?;

    return Ok(());
}
