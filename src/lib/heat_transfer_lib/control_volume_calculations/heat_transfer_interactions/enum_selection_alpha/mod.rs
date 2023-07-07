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
use std::f64::consts::PI;

use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::CVType::*;


use crate::heat_transfer_lib::control_volume_calculations
::heat_transfer_interactions::*;

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
            caclulate_between_two_singular_cv_nodes(
                single_cv_1, 
                single_cv_2, 
                interaction),
        (ArrayCV, SingleCV(_)) =>
            Err("Array CV calcs not yet implemented".to_string()),
        (SingleCV(_), ArrayCV) =>
            Err("Array CV calcs not yet implemented".to_string()),
        (ArrayCV, ArrayCV) =>
            Err("Array CV calcs not yet implemented".to_string()),
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
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_control_volume_boundary_condition_serial(
    control_vol: &mut CVType,
    boundary_condition: &mut BCType,
    interaction: HeatTransferInteractionType) -> Result<(),String> {

    let cv_bc_result = match (control_vol, boundary_condition) {
        (SingleCV(single_cv), BCType::UserSpecifiedHeatAddition(heat_rate)) 
            => calculate_single_cv_node_constant_heat_addition(
                single_cv, *heat_rate, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedHeatFlux(heat_flux))
            => calculate_single_cv_node_constant_heat_flux(
                single_cv, *heat_flux, interaction),
        (SingleCV(single_cv), BCType::UserSpecifiedTemperature(bc_temperature))
            => calculate_single_cv_node_constant_temperature(
                single_cv, *bc_temperature, interaction),
        (ArrayCV,_) => Err("array cv not yet implemented".to_string()),
    };

    return cv_bc_result;
}


/// suppose the control volume interacts with a BC which is 
/// a constant heat addition, which is a constant power rating
/// 
///
/// cooling is also possible, just supply a negative Power 
/// quantity
///
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_single_cv_node_constant_heat_addition(
    control_vol: &mut SingleCVNode,
    heat_added_to_control_vol: Power,
    interaction: HeatTransferInteractionType
    ) -> Result<(), String> {

    // ensure that the interaction is UserSpecifiedHeatAddition
    // otherwise, return error 

    match interaction {
        heat_transfer_lib::control_volume_calculations::
            heat_transfer_interactions::
            enums_alpha::HeatTransferInteractionType::
            UserSpecifiedHeatAddition => {
                // return a void value, that would be dropped 
                // instantly
                //
                // it pretty much has the same meaning as break
                ()
            },

        _ => return Err("you need to specify that the interaction type \n 
            is UserSpecifiedHeatAddition".to_string()),
    };


    control_vol.rate_enthalpy_change_vector.
        push(heat_added_to_control_vol);

    return Ok(());
}

/// suppose control volume interacts with a constant heat flux BC
/// we will need a sort of surface area in order to determine the 
/// power added
fn calculate_single_cv_node_constant_heat_flux(
    control_vol: &mut SingleCVNode,
    heat_flux_into_control_vol: HeatFluxDensity,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

     let heat_transfer_area: Area = match interaction{
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            DualCartesianThermalConductance(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(_, _, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or \n
                Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

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

    };

    let heat_flowrate_into_control_vol: Power = 
    heat_flux_into_control_vol * heat_transfer_area;

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_into_control_vol);

    return Ok(());


}

/// suppose control volume interacts with constant temperature BC
/// we need some thermal conductance value to obtain a power
/// value
fn calculate_single_cv_node_constant_temperature(
    control_vol: &mut SingleCVNode,
    boundary_condition_temperature: ThermodynamicTemperature,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    // first let's get the control volume temperatures out
    
    let cv_enthalpy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    let cv_material = control_vol.material_control_volume;

    let cv_pressure = control_vol.pressure_control_volume;

    let cv_temperature = heat_transfer_lib::thermophysical_properties:: 
        specific_enthalpy::temperature_from_specific_enthalpy(
            cv_material, 
            cv_enthalpy, 
            cv_pressure)?;

    // for now we assume the boundary condition pressure is the same 
    // as the control volume pressure, because pressure is not 
    // specified or anything

    let bc_pressure = cv_pressure.clone();
    
    // we'll need thermal conductance 

    let cv_bc_conductance: ThermalConductance = 
    get_thermal_conductance(
        cv_temperature, 
        boundary_condition_temperature, 
        cv_pressure, 
        bc_pressure, 
        interaction)?;

    // with conductance settled, we should be able to calculate a power 

    let bc_temp_minus_cv_temp_kelvin: f64 = 
        boundary_condition_temperature.get::<kelvin>() - 
        cv_temperature.get::<kelvin>();

    let bc_temp_mins_cv_temp: TemperatureInterval = 
    TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(
        bc_temp_minus_cv_temp_kelvin);

    // heat flow to the destination 
    // be is proportional to -(T_final - T_initial) 
    // or -(T_destination - T_source)
    //
    // so if the destination is the boundary condition,
    //
    let heat_flowrate_from_cv_to_bc: Power = 
    - cv_bc_conductance * bc_temp_mins_cv_temp;

    // now, push the power change or heat flowrate 
    // to the control volume 
    //

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_from_cv_to_bc);

    // and we done!
    return Ok(());
}

// this function specifically calculates interaction 
// between two single CV nodes
// 
#[inline]
fn caclulate_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType)-> Result<(), String>{

    // let's get the two temperatures of the control volumes first
    // so let me get the enthalpies, and then their respective 
    // temperatures 

    let single_cv_1_enthalpy = single_cv_1.
        current_timestep_control_volume_specific_enthalpy.clone();
    let single_cv_2_enthalpy = single_cv_2.
        current_timestep_control_volume_specific_enthalpy.clone();

    // to get the temperatures, we'll need the material as well 
    let single_cv_1_material = single_cv_1.material_control_volume.clone();
    let single_cv_2_material = single_cv_2.material_control_volume.clone();

    // we'll also need to get their pressures 
    let single_cv_1_pressure = single_cv_1.pressure_control_volume.clone();
    let single_cv_2_pressure = single_cv_2.pressure_control_volume.clone();

    // we will now get their respective temperatures 
    let single_cv_1_temperature = temperature_from_specific_enthalpy(
        single_cv_1_material, 
        single_cv_1_enthalpy, 
        single_cv_1_pressure)?;
    let single_cv_2_temperature = temperature_from_specific_enthalpy(
        single_cv_2_material, 
        single_cv_2_enthalpy, 
        single_cv_2_pressure)?;

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
    // // TODO: probably change the unwrap for later
    let thermal_conductance = get_thermal_conductance(
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

    return Ok(());

}

/// this is thermal conductance function 
/// it calculates thermal conductance based on the supplied enum
///
/// TODO: probably want to test this function out
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
                }
                ,
            _ => return Err("define other enum options".to_string()),
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
