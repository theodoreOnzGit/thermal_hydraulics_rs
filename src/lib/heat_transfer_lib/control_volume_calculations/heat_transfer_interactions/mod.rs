use std::f64::consts::PI;

use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{f64::*, pressure::atmosphere};
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use crate::heat_transfer_lib::thermophysical_properties
::{Material, specific_enthalpy::specific_enthalpy};
use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::HeatTransferEntity;

/// This is a calculation library for all interactions between 
/// the various heat transfer entities
pub(in crate::heat_transfer_lib::
    control_volume_calculations::heat_transfer_interactions) 
mod calculations;
use calculations::*;

use super::common_functions::obtain_power_two_convection_one_conduction_thermal_resistance;

/// basically an enum for you to specify 
/// if the liquid on the inner curved surface of the shell or outer 
/// curved surface of the shell
///
/// in the context of a convection and conductivity 
/// thermal resistance calculation,
///
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum CylindricalAndSphericalSolidFluidArrangement {
    /// indicates that fluid in the inner side of a curved shell 
    ///
    /// -----------------------------------------> r
    /// fluid               ||                  solid
    ///
    FluidOnInnerSurfaceOfSolidShell,
    /// indicates that fluid in the inner side of a curved shell 
    ///
    /// -----------------------------------------> r
    /// solid               ||                  fluid
    ///
    FluidOnOuterSurfaceOfSolidShell
}



/// Contains possible heat transfer interactions between the nodes
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum HeatTransferInteractionType {
    /// The user specifies a thermal conductance between the nodes
    /// in units of power/kelvin
    UserSpecifiedThermalConductance(ThermalConductance),

    /// 1D Cartesian Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have two control volumes, each node represents a control 
    /// volume
    ///
    /// // ----------------------------
    /// // |                          |
    /// // *                          *
    /// // |                          |
    /// // ----------------------------
    /// // cv_1                      cv_2
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have one material which determines conductivity 
    /// and then a length which determines the distance between 
    /// the two control volumes
    ///
    SingleCartesianThermalConductanceOneDimension(
        Material,
        XThicknessThermalConduction
    ),

    /// 1D Cartesian Coordinates Thermal Resistance, for solids only
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // cv_1                      cv_2                     cv_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two control volumes
    ///
    /// Information must be passed in as a tuple,
    ///
    ///
    DualCartesianThermalConductance(
        (XThicknessThermalConduction, Length),
        (XThicknessThermalConduction, Length),
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes 
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // cv_1                      cv_2                     cv_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two control volumes 
    ///
    /// one also needs to determine the 
    /// inner diameter, outer diameter and length of the tube 
    /// 
    ///
    /// 
    ///
    DualCylindricalThermalConductance(
        (Material,RadialCylindricalThicknessThermalConduction),
        (Material,RadialCylindricalThicknessThermalConduction),
        (InnerDiameterThermalConduction, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction
    )
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes along the outer wall
    ///
    /// -------------------------------------------------------> r
    /// // ----------------------------
    /// // |                          |                          
    /// // * solid_cv_1               *                          *
    /// // |                          |                         (T_f) 
    /// // ----------------------------
    /// //                        solid_surface              Fluid_node
    ///
    /// Where r is the radius 
    /// basically the liquid is on the outside (larger r)
    ///
    /// between solid_cv_1 and the solid_surface 
    /// cv_2 there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// between solid_surface and fluid_node, there is convection resistance
    /// specified by a Nusselt Number so that we get a heat transfer 
    /// coefficient
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two control volumes
    ///
    /// the thermal conductance is determined by 
    /// Thermal conductance 
    /// /// (2 * pi * L * K)/
    /// ln(outer_radius/inner_radius)
    ///
    /// using obtain_thermal_conductance_annular_cylinder
    /// under common_functions
    ///
    ///
    /// For convection, the heat flux from solid surface to fluid 
    /// is:
    ///
    /// q = h A(T_s - T_f)
    /// 
    /// for hA
    /// surface area is calculated by specifying an outer diameter 
    /// and a cylindrical axial length
    ///
    ///
    ///
    CylindricalConductionConvectionLiquidOutside(
        (Material,RadialCylindricalThicknessThermalConduction,
        ThermodynamicTemperature,Pressure),
        (HeatTransfer, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction),
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes along the outer wall
    ///
    /// -------------------------------------------------------> r
    /// //                           ----------------------------
    /// //                           |                          |                          
    /// // *                         *         solid_cv_1       *                  
    /// //                           |                          |                   
    /// // fluid node                ----------------------------
    /// // (T_f)                solid_surface
    ///
    /// Where r is the radius 
    /// basically the liquid is on the inside (larger smaller r)
    ///
    /// between solid_cv_1 and solid_surface 
    /// there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// between solid_surface and fluid_node, there is convection resistance
    /// specified by a Nusselt Number so that we get a heat transfer 
    /// coefficient
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two control volumes
    ///
    /// the thermal conductance is determined by 
    /// Thermal conductance 
    /// /// (2 * pi * L * K)/
    /// ln(outer_radius/inner_radius)
    ///
    /// using obtain_thermal_conductance_annular_cylinder
    /// under common_functions
    ///
    ///
    /// For convection, the heat flux from solid surface to fluid 
    /// is:
    ///
    /// q = h A(T_s - T_f)
    /// 
    /// for hA
    /// surface area is calculated by specifying an outer diameter 
    /// and a cylindrical axial length
    ///
    ///
    ///
    CylindricalConductionConvectionLiquidInside(
        (Material,RadialCylindricalThicknessThermalConduction,
        ThermodynamicTemperature,Pressure),
        (HeatTransfer, 
        InnerDiameterThermalConduction, 
        CylinderLengthThermalConduction),
    ),


    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}

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
pub fn link_heat_transfer_entity(entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType){

    // first thing first, probably want to unpack the enums to obtain 
    // the underlying control volume and BCs

    let control_vol_1: Option<&mut CVType> = match entity_1 {
        HeatTransferEntity::ControlVolume(control_vol_type) 
            => Some(control_vol_type),
        _ => None,
    };

    let control_vol_2: Option<&mut CVType> = match entity_2 {
        HeatTransferEntity::ControlVolume(control_vol_type) 
            => Some(control_vol_type),
        _ => None,
    };

    // I'll pass in both these option types into a function which 
    // calculates specifically for two control volumes


   
}

// the job of this function is to take in a control volume 
// and then mutate it by calculating its interaction
#[inline]
fn calculate_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType){

    // let me first match my control volumes to their various types
    //
    // // TODO: handle the case of heat addition 
    //
    // // TODO: handle the case of conduction

}

// this function specifically interacts between two single CV nodes
//
#[inline]
fn caclulate_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType){

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
        single_cv_1_pressure)
        .unwrap();
    let single_cv_2_temperature = temperature_from_specific_enthalpy(
        single_cv_2_material, 
        single_cv_2_enthalpy, 
        single_cv_2_pressure)
        .unwrap();

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
        interaction).unwrap();

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
    single_cv_1.rate_enthalpy_change_vector.
        push(heat_flowrate_from_cv_1_to_cv_2);

}

/// this is thermal conductance function 
/// it calculates thermal conductance based on the supplied enum
///
/// TODO: probably want to test this function out
fn get_thermal_conductance(
    material_temperature_1: ThermodynamicTemperature,
    material_temperature_2: ThermodynamicTemperature,
    material_pressure_1: Pressure,
    material_pressure_2: Pressure,
    interaction: HeatTransferInteractionType) -> Result<ThermalConductance, String> 
{

    let conductance: ThermalConductance = match 
        interaction {
            HeatTransferInteractionType::UserSpecifiedThermalConductance(
                user_specified_conductance) => user_specified_conductance,
            HeatTransferInteractionType
                ::SingleCartesianThermalConductanceOneDimension(
                material,thickness) => get_conductance_single_cartesian_one_dimension(
                    material,
                    material_temperature_1, 
                    material_temperature_2, 
                    material_pressure_1, 
                    material_pressure_2, 
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
            ) => get_conductance_cylindrical_radial_two_materials(
                    inner_material,
                    outer_material,
                    material_temperature_1, //convention, 1 is inner shell
                    material_temperature_2, // convention 2, is outer shell
                    material_pressure_1,
                    material_pressure_2,
                    inner_diameter,
                    inner_shell_thickness,
                    outer_shell_thickness,
                    cylinder_length,
                )?,
            _ => todo!("define other enum options"),
        };

    return Ok(conductance);
}

/// todo, test to be implemented
#[test]
fn thermal_conductance_test_convection_conduction_boundary(){

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




///
/// For this part, we determine how heat_transfer_entities interact 
/// with each other by using a function 
/// The function will take in a HeatTransferInteractionType enum
/// which you must first initiate
///
/// Notes for parallel computation: 
///
/// if we use mutable borrows early, it will be difficult to use 
/// this in a multithreaded manner
///
/// We can use interior mutability to allow for mutex locks and Arc 
/// types to be used so that parallel calculations can be performed 
///
/// lastly, we can just make sure we return a value from the placeholder 
/// function, and then manually map the heat transfer value to each 
/// heat_transfer_entities pair 
fn placeholder_function_suitable_for_parallel_computation(){}
