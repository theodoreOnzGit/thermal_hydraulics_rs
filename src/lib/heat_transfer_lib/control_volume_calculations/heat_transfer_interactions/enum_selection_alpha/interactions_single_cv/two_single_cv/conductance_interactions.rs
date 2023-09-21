use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::try_get_temperature_from_h;



use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
#[inline]
pub fn calculate_conductance_interaction_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType)-> Result<(), String>{

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

    let experimental_code = false; 
    if experimental_code {

        single_cv_1_temperature = try_get_temperature_from_h(
            single_cv_1_material, 
            single_cv_1_enthalpy, 
            single_cv_1_pressure)?;
        single_cv_2_temperature = try_get_temperature_from_h(
            single_cv_2_material, 
            single_cv_2_enthalpy, 
            single_cv_2_pressure)?;

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

