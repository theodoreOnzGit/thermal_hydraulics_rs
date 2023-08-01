
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

use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
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
        Solid(_) => {
            ()
        },
        Liquid(_) => {
            todo!("need to calculate convection based time scales")
        },
    }

    return Ok(());
}

/// suppose control volume interacts with a constant heat flux BC
/// we will need a sort of surface area in order to determine the 
/// power added
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_single_cv_node_constant_heat_flux(
    control_vol: &mut SingleCVNode,
    heat_flux_into_control_vol: HeatFluxDensity,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    // first, obtain a heat transfer area from the constant heat flux 
    // BC
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

        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(_) => 
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

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
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
        Solid(_) => {
            ()
        },
        Liquid(_) => {
            todo!("need to calculate convection based time scales")
        },
    }
    
    return Ok(());


}

/// suppose control volume interacts with constant temperature BC
/// we need some thermal conductance value to obtain a power
/// value
pub (in crate::heat_transfer_lib::control_volume_calculations)
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
        push(-heat_flowrate_from_cv_to_bc);

    // for constant temperature BC in interaction with CV,
    // we only need take into consideration the CV 
    // timescale, this is based on the fourier number or Nusselt 
    // number edited fourier number
    //
    // for solids mesh fourier number need only 
    // be done once, not every time 
    // an interaction is formed 
    //
    // probably the cell stability fourier number will be done in the 
    // constructor. however, with convection, the time scale must be 
    // recalculated at every time step. so it really depends whether 
    // it's solid or fluid control volume
    // 

    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    match cv_material {
        Solid(_) => {
            ()
        },
        Liquid(_) => {
            todo!("need to calculate convection based time scales")
        },
    }
    return Ok(());
}

// this function specifically calculates interaction 
// between two single CV nodes
// 
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations)
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
    match single_cv_1_material {
        Solid(_) => {
            ()
        },
        Liquid(_) => {
            todo!("need to calculate convection based time scales")
        },
    }

    match single_cv_2_material {
        Solid(_) => {
            ()
        },
        Liquid(_) => {
            todo!("need to calculate convection based time scales")
        },
    }
    return Ok(());

}
