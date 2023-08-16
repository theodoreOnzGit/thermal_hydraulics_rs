//! enum selection basically handles how control volumes (CVs) interact 
//! with each other and with boundary conditions (BCs)
//! 
//! together, the CVs and BCs are part of a superset or enum called 
//! heat_transfer_entities
//!
//! This enum selection decides what happens if two CVs interact with 
//! each other, or a CV and BC interact with each other, or if two BCs 
//! interact with each other 
//!
//! CVs differ from BCs in that BCs do not have thermal inertia 
//! or any heat capacity. BCs will supply a heat or stay at a fixed 
//! Temperature (thermal reservoir)
//!
//!
//! Note that if two BCs interact with each other, it kind of doesn't 
//! matter because you cannot change temperatures of each BC
//!



use uom::si::f64::*;

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::CVType::*;


use crate::heat_transfer_lib::control_volume_calculations
::heat_transfer_interactions::*;

pub (in crate::heat_transfer_lib::control_volume_calculations) 
    mod interactions_single_cv;
use interactions_single_cv::*;

pub (in crate::heat_transfer_lib::control_volume_calculations) 
    mod timestep_control_single_cv;
use timestep_control_single_cv::*;

// the job of this function is to take in a control volume 
// and then mutate it by calculating its interaction
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<(),String>{

    // let me first match my control volumes to their various types
    // at each matching arm, use those as inputs

    let cv_result = match (control_vol_1, control_vol_2) {
        (SingleCV(single_cv_1), SingleCV(single_cv_2)) =>
            calculate_between_two_singular_cv_nodes(
                single_cv_1, 
                single_cv_2, 
                interaction),
        // basically like this 
        //
        // (array_cv_back -----array_cv --- array_cv_front) --- (single_cv) 
        (ArrayCV(array_cv), SingleCV(single_cv)) =>
            {
                array_cv.link_single_cv_to_higher_side(
                    single_cv,
                    interaction)
            },
        // basically like this 
        //
        // (single_cv) --- (array_cv_back -----array_cv --- array_cv_front)
        (SingleCV(single_cv), ArrayCV(array_cv)) =>
            {
                array_cv.link_single_cv_to_lower_side(
                    single_cv,
                    interaction)
            },
        // basically like this 
        //
        // (back --- cv_1 --- front) ---- (back --- cv_2 --- front)
        (ArrayCV(array_cv_1), ArrayCV(array_cv_2)) =>
            {
                array_cv_1.link_array_cv_to_the_front_of_this_cv(
                    array_cv_2,
                    interaction)
            },
    };


    return cv_result;

}

/// the job of this function is to calculate time step 
/// when two control volumes interact, the control volumes can be 
/// single control volumes or array control volumes
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_timestep_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<Time,String>{

    // let me first match my control volumes to their various types
    // at each matching arm, use those as inputs

    let cv_result = match (control_vol_1, control_vol_2) {
        (SingleCV(single_cv_1), SingleCV(single_cv_2)) =>
            calculate_mesh_stability_timestep_for_two_single_cv_nodes(
                single_cv_1, 
                single_cv_2, 
                interaction),
        (ArrayCV(array_cv_type), SingleCV(single_cv_node)) =>
            {
                array_cv_type.
                    calculate_timestep_for_single_cv_to_front_of_array_cv(
                        single_cv_node,
                        interaction)
            },
        (SingleCV(single_cv_node), ArrayCV(array_cv_type)) =>
            {
                array_cv_type.
                    calculate_timestep_for_single_cv_to_back_of_array_cv(
                        single_cv_node,
                        interaction)
            },
        (ArrayCV(array_cv_type_self), ArrayCV(array_cv_other)) =>
            {
                array_cv_type_self.
                    calculate_timestep_for_array_cv_to_front_of_this_array_cv(
                        array_cv_other,
                        interaction)
            },
    };


    return cv_result;

}

/// the job of this function is to handle interactions between 
/// boundary conditions 
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_boundary_condition_serial(
    _boundary_condition_1: &mut BCType,
    _boundary_condition_2: &mut BCType,
    _interaction: HeatTransferInteractionType) -> Result<(),String>{


    return Err("interactions between two boundary conditions \n
        not implemented".to_string());
}

/// the job of this function is to calcualte timestep between 
/// boundary conditions 
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_timestep_boundary_condition_serial(
    _boundary_condition_1: &mut BCType,
    _boundary_condition_2: &mut BCType,
    _interaction: HeatTransferInteractionType) -> Result<Time,String>{


    return Err("timesteps between two boundary conditions \n
        not implemented".to_string());
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
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn attach_boundary_condition_to_control_volume_front_serial(
    control_vol: &mut CVType,
    boundary_condition: &mut BCType,
    interaction: HeatTransferInteractionType) -> Result<(),String> {

    let cv_bc_result = match (control_vol, boundary_condition) {
        (SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(heat_rate)) 
            => calculate_constant_heat_addition_front_single_cv_back(
                single_cv, *heat_rate, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(heat_flux))
            => calculate_constant_heat_flux_front_single_cv_back(
                single_cv, *heat_flux, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedTemperature(bc_temperature))
            => calculate_constant_temperature_front_single_cv_back(
                single_cv, *bc_temperature, interaction),
        (ArrayCV(cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            cv.link_heat_flux_bc_to_front_of_this_cv(
                *heat_flux,
                interaction)
        },
        (ArrayCV(cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            cv.link_heat_addition_to_front_of_this_cv(
                *heat_rate,
                interaction)
        },
        (ArrayCV(cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            cv.link_constant_temperature_to_front_of_this_cv(
                *bc_temperature,
                interaction)
        },
    };

    return cv_bc_result;
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
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn attach_boundary_condition_to_control_volume_back_serial(
    boundary_condition: &mut BCType,
    control_vol: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<(),String> {

    let cv_bc_result = match (control_vol, boundary_condition) {
        (SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(heat_rate)) 
            => calculate_single_cv_front_constant_heat_addition_back(
                *heat_rate, single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(heat_flux))
            => calculate_single_cv_front_heat_flux_back(
                *heat_flux,single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedTemperature(bc_temperature))
            => calculate_single_cv_node_front_constant_temperature_back(
                *bc_temperature, single_cv, interaction),
        (ArrayCV(cv),BCType::UserSpecifiedHeatFlux(heat_flux)) => {
            cv.link_heat_flux_bc_to_back_of_this_cv(
                *heat_flux,
                interaction)

        },
        (ArrayCV(cv),BCType::UserSpecifiedHeatAddition(heat_rate)) => {
            cv.link_heat_addition_to_back_of_this_cv(
                *heat_rate,
                interaction)
        },
        (ArrayCV(cv),BCType::UserSpecifiedTemperature(bc_temperature)) => {
            cv.link_constant_temperature_to_back_of_this_cv(
                *bc_temperature,
                interaction)
        },
    };

    return cv_bc_result;
}


/// calculates timestep for control voluem and boundary condition 
/// in serial
/// for arrayCVs, the boundary condition is attached to the front 
/// of the control volume
///
/// (back --- cv_1 --- front) ---- (boundary condition)
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_timestep_control_volume_front_to_boundary_condition_serial(
    control_vol: &mut CVType,
    boundary_condition: &mut BCType,
    interaction: HeatTransferInteractionType) -> Result<Time,String> {

    // to prevent myself from too much boiler plate code,
    // i'm nesting a closure here 
    //
    // this will extract the front control volume (or outer control 
    // volume) and calculate the timestep based on that
    let match_array_cv_front = |array_cv_type: &mut ArrayCVType|{

        match array_cv_type {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                // get a mutable reference for the front cv 

                let front_cv_reference = 
                &mut cartesian_array_cv.outer_single_cv;


                // once that is done, use the same function
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                    front_cv_reference, interaction)

            },
        }
    };

    let cv_bc_result = match (control_vol, boundary_condition) {
        (SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(_heat_rate)) 
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(_heat_flux))
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedTemperature(_bc_temperature))
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (ArrayCV(array_cv_type),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 
            match_array_cv_front(array_cv_type)
        },
        (ArrayCV(array_cv_type),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 

            match_array_cv_front(array_cv_type)
        },
        (ArrayCV(array_cv_type),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            // match the cv to the cv type, extract the back boundary 
            // condition 

            match_array_cv_front(array_cv_type)
        },
    };

    return cv_bc_result;
}

/// calculates timestep for control voluem and boundary condition 
/// in serial
/// for arrayCVs, the boundary condition is attached to the back
/// of the control volume
///
/// (boundary condition) ----- (back --- cv_1 --- front) 
///
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_timestep_control_volume_back_to_boundary_condition_serial(
    boundary_condition: &mut BCType,
    control_vol: &mut CVType,
    interaction: HeatTransferInteractionType) -> Result<Time,String> {

    // to prevent myself from too much boiler plate code,
    // i'm nesting a closure here 
    //
    // this will extract the back control volume (or outer control 
    // volume) and calculate the timestep based on that
    let match_array_cv_back = |array_cv_type: &mut ArrayCVType|{

        // match the cv to the cv type, extract the back boundary 
        // condition 

        match array_cv_type {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                // get a mutable reference for the back cv 

                let back_cv_reference = 
                &mut cartesian_array_cv.inner_single_cv;


                // once that is done, recycle the condition
                calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                    back_cv_reference, interaction)

            },
        }
    };

    let cv_bc_result = match (control_vol, boundary_condition) {
        (SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(_heat_rate)) 
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(_heat_flux))
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedTemperature(_bc_temperature))
            => calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                single_cv, interaction),
        (ArrayCV(array_cv_type),BCType::UserSpecifiedHeatFlux(_heat_flux)) => {
            match_array_cv_back(array_cv_type)
        },
        (ArrayCV(array_cv_type),BCType::UserSpecifiedHeatAddition(_heat_rate)) => {
            match_array_cv_back(array_cv_type)
        },
        (ArrayCV(array_cv_type),BCType::UserSpecifiedTemperature(_bc_temperature)) => {
            match_array_cv_back(array_cv_type)
        },
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
fn get_thermal_conductance(
    temperature_1: ThermodynamicTemperature,
    temperature_2: ThermodynamicTemperature,
    pressure_1: Pressure,
    pressure_2: Pressure,
    interaction: HeatTransferInteractionType) 
-> Result<ThermalConductance, String> 
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


                        return Err(error_str
                        );
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
                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());
                },

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCustomArea(_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(_,_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalInnerArea(_,_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

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
                return Err("advection interaction types \n 
                do not correspond to conductance".to_string());
            },

        };

    return Ok(conductance);
}

/// basically have a simple test here, I want to  check whether 
/// the thermal conductance returned by the enum system works 
///
#[test]
fn thermal_conductance_test_convection_conduction_boundary(){

    use uom::si::{f64::*, pressure::atmosphere};
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use crate::heat_transfer_lib::thermophysical_properties::
    thermal_conductivity::thermal_conductivity;
    use crate::heat_transfer_lib::thermophysical_properties::Material;
    use std::f64::consts::PI;
    // first we can use the existing common_functions 
    // model to obtain power 
    // between a one convection one conduction thermal resistance
    // let's suppose we have a pipe of 0.015 m inner diameter 
    // 0.020 m outer diameter 
    // this is a cylindrical surface
    // pipe is 1m long
    //
    // pipe is made of copper
    //
    // the pipe itself has a heat transfer coefficient of 
    // 100 W/(m^2 K) inside and 20 W/(m^2 K) 
    //
    // The ambient air temperature is 20 C, whereas the 
    // fluid temperature is 100C 

    let temperature_of_heat_recipient = ThermodynamicTemperature::new
        ::<uom::si::thermodynamic_temperature::degree_celsius>(20.0);

    let temperature_of_heat_source = ThermodynamicTemperature::new
        ::<uom::si::thermodynamic_temperature::degree_celsius>(100.0);

    // now we define the inner surface,
    // area is pi* id * L (inner diameter)
    // this is boundary between pipe and liquid

    let heat_transfer_coefficient_1 = HeatTransfer::new 
        ::<watt_per_square_meter_kelvin>(100.0);

    let id = Length::new::<uom::si::length::meter>(0.015);
    let pipe_length = Length::new::<uom::si::length::meter>(1.0);

    let average_surface_area_1: Area = id * pipe_length * PI;

    // then define the outer surface, this is boundary between pipe and 
    // ambient air

    let heat_transfer_coefficient_2 = HeatTransfer::new 
    ::<watt_per_square_meter_kelvin>(20.0);

    let od = Length::new::<uom::si::length::meter>(0.020);
    let pipe_length = Length::new::<uom::si::length::meter>(1.0);

    let average_surface_area_2: Area = od * pipe_length * PI;

    // now for the conductance thermal resistance, while there is 
    // a formal integral with which to get average surface area for
    // the thermal conductance, perhaps with some use of log mean area 
    // or something like that, I'm just going to lazily average the 
    // surface area to get a ballpark figure 
    
    let lazily_averaged_surface_area_for_conduction: Area = 
    0.5*(average_surface_area_1 + average_surface_area_2);

    // now let's get copper thermal conductivity at some mean temperature,
    // say 70C, of course, one can calculate this out iteratively 
    // but i only ballpark 

    let copper_temperature = ThermodynamicTemperature::new
        ::<uom::si::thermodynamic_temperature::degree_celsius>(70.0);

    let copper: Material = Material::Solid(
        crate::heat_transfer_lib::
        thermophysical_properties::SolidMaterial::Copper);

    // now let's obtain thermal conductivity

    let copper_pressure = Pressure::new::<atmosphere>(1.0);

    let average_thermal_conductivity = 
    thermal_conductivity(copper, copper_temperature, copper_pressure)
        .unwrap();

    // now let's obtain the power, 
    use crate::heat_transfer_lib::control_volume_calculations::common_functions::
    obtain_power_two_convection_one_conduction_thermal_resistance;
    let ballpark_estimate_power: Power = 
    obtain_power_two_convection_one_conduction_thermal_resistance(
        temperature_of_heat_recipient, 
        temperature_of_heat_source, 
        average_surface_area_1, 
        heat_transfer_coefficient_1, 
        average_surface_area_2, 
        heat_transfer_coefficient_2, 
        average_thermal_conductivity, 
        lazily_averaged_surface_area_for_conduction, 
        pipe_length);

    // the power across this is about 75 watts
    approx::assert_relative_eq!(
        ballpark_estimate_power.value,
        75.6831308452296,
        epsilon=0.01);

    // the conductance (Htc) is about 
    // q = Htc (T_1-T_2)

    let temperature_diff: TemperatureInterval = TemperatureInterval::new
        ::<uom::si::temperature_interval::degree_celsius>(80.0);

    let thermal_conductance_ballpark_value: ThermalConductance = 
    ballpark_estimate_power/temperature_diff;

    // conductance is about 0.946 watt/kelvin
    approx::assert_relative_eq!(
        thermal_conductance_ballpark_value.value,
        0.94603913556537,
        epsilon=0.01);

    // now then, let's get conductance from here.

    // first we need an interaction type: 
    // and two material temperatures and pressures 
    // we assume copper is 70C for simplicity

    // we'll need to do two halves, one is the copper 
    // with hot liquid inside 
    // the help messages aren't too helpful here
    // from the compiler
    // may want to insert a struct instead

    // the thickness I expect here is perhaps 1/4 of the pipe thickness 
    // this is a rough estimate

    let radial_thicnkess = (od-id)/4.0;
    let radial_thicnkess: RadialCylindricalThicknessThermalConduction = 
    RadialCylindricalThicknessThermalConduction::from(radial_thicnkess);


    // it's quite a clunky way to define geometry and thermal resistance 
    // but it will have to do for now
    //
    // moreover, one has to decide where best to situate the control 
    // vol nodes so that thermal resistance is best defined
    let copper_hot_liquid_interaction = 
        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
                (// kind of clunky to use but ok
                    copper,
                    radial_thicnkess,
                    copper_temperature,
                    copper_pressure,
                ), 
                (
                    heat_transfer_coefficient_1,
                    InnerDiameterThermalConduction::from(id),
                    CylinderLengthThermalConduction::from(pipe_length),
                )
            );

    let inner_convection_boundary_conductance: ThermalConductance = 
    get_thermal_conductance(
        copper_temperature, 
        copper_temperature, 
        copper_pressure, //these don't really need to be defined, 
        //probably fine tune later
        copper_pressure, 
        copper_hot_liquid_interaction).unwrap();

    let inner_convection_boundary_thermal_resistance = 
    1.0/inner_convection_boundary_conductance;

    // now let's do it for the outer boundary 

    let copper_air_interaction = HeatTransferInteractionType::
        CylindricalConductionConvectionLiquidOutside(
            (
                copper,
                radial_thicnkess,
                copper_temperature,
                copper_pressure,
            ), 
            (
                heat_transfer_coefficient_2,
                OuterDiameterThermalConduction::from(od),
                CylinderLengthThermalConduction::from(pipe_length),
            )
        );

    let outer_convection_boundary_conductance: ThermalConductance = 
    get_thermal_conductance(
        copper_temperature, 
        copper_temperature, 
        copper_pressure, //these don't really need to be defined, 
        //probably fine tune later
        copper_pressure, 
        copper_air_interaction).unwrap();

    let outer_convection_boundary_thermal_resistance = 
    1.0/outer_convection_boundary_conductance;

    // now let's do it for the middle boundary 
    // where there are two layers of copper
    
    let inner_interim_diameter: Length = id + radial_thicnkess.into();
    let outer_interim_diameter: Length = od - radial_thicnkess.into();

    let copper_shell_internal_interaction = 
    HeatTransferInteractionType::
        DualCylindricalThermalConductance(
            (copper, radial_thicnkess), 
            (copper, radial_thicnkess), 
            (
                InnerDiameterThermalConduction::from(inner_interim_diameter),
                OuterDiameterThermalConduction::from(outer_interim_diameter),
                CylinderLengthThermalConduction::from(pipe_length),
            )
        );

    let copper_shell_inside_thermal_conductance: ThermalConductance = 
    get_thermal_conductance(
        copper_temperature, 
        copper_temperature, 
        copper_pressure, 
        copper_pressure, 
        copper_shell_internal_interaction
    ).unwrap();

    let copper_shell_inside_thermal_resistance = 
    1.0/copper_shell_inside_thermal_conductance;

    // now find overall conductance 

    let overall_resistance = copper_shell_inside_thermal_resistance +
    outer_convection_boundary_thermal_resistance +
    inner_convection_boundary_thermal_resistance;

    let overall_conductance: ThermalConductance = 
    1.0/overall_resistance;

    // values match to within 0.5%
    // yes!
    //
    // must say though, it's kind of user unfriendly to use
    // for now, i'll leave it, it seems to work
    approx::assert_relative_eq!(
        thermal_conductance_ballpark_value.value,
        overall_conductance.value,
        epsilon=0.005);

}
