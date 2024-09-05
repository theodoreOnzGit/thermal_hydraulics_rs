
#[cfg(test)]
/// function to verify the dhx analytical solution
/// version 1,
/// the DHX here uses uncalibrated Gnielinski correlations 
/// to estimate heat transfer coefficients
pub fn verify_coupled_dhx_analytical_solution_version_1(
    input_power_watts: f64,
    tchx_outlet_temperature_set_point_degc: f64,
    experimental_dracs_mass_flowrate_kg_per_s: f64,
    experimental_primary_mass_flowrate_kg_per_s: f64) -> 
Result<(),crate::thermal_hydraulics_error::ThermalHydraulicsLibError>{
    use uom::si::{f64::*, mass_rate::kilogram_per_second, power::watt};

    use uom::si::{frequency::hertz, ratio::ratio, time::millisecond};

    use crate::boussinesq_solver::pre_built_components::ciet_isothermal_test_components::*;
    use crate::boussinesq_solver::pre_built_components::ciet_steady_state_natural_circulation_test_components::coupled_dracs_loop_tests::dhx_constructor::new_dhx_sthe_version_1;
    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::dracs_loop_components::*;
    use uom::ConstZero;

    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use uom::si::time::second;

    let input_power = Power::new::<watt>(input_power_watts);
    let _experimental_dracs_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_dracs_mass_flowrate_kg_per_s);
    let _experimental_primary_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_primary_mass_flowrate_kg_per_s);

    let tchx_outlet_temperature_set_point = 
        ThermodynamicTemperature::new::<degree_celsius>(
            tchx_outlet_temperature_set_point_degc);
    use chem_eng_real_time_process_control_simulator::alpha_nightly::transfer_fn_wrapper_and_enums::TransferFnTraits;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::ProportionalController;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::AnalogController;

    // max error is 0.5% according to SAM 
    // is okay, because typical flowmeter measurement error is 2% anyway
    let timestep = Time::new::<second>(0.5);
    let heat_rate_through_dhx = input_power;
    let mut tchx_heat_transfer_coeff: HeatTransfer;

    let reference_tchx_htc = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(40.0);
    let average_temperature_for_density_calcs = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);
    // let's calculate 3800 seconds of simulated time 
    // it takes about that long for the temperature to settle down
    // this is compared to value at 4000s

    let mut current_simulation_time = Time::ZERO;
    let max_simulation_time = Time::new::<second>(3800.0);

    // PID controller settings
    let controller_gain = Ratio::new::<ratio>(1.75);
    let integral_time: Time = controller_gain / Frequency::new::<hertz>(1.0);
    let derivative_time: Time = Time::new::<second>(1.0);
    // derivative time ratio
    let alpha: Ratio = Ratio::new::<ratio>(1.0);

    let mut pid_controller: AnalogController = 
        AnalogController::new_filtered_pid_controller(controller_gain,
            integral_time,
            derivative_time,
            alpha).unwrap();

    // we also have a measurement delay of 0.0001 s 
    // or 0.1 ms
    let measurement_delay = Time::new::<millisecond>(0.1);

    let mut measurement_delay_block: AnalogController = 
        ProportionalController::new(Ratio::new::<ratio>(1.0)).unwrap().into();

    measurement_delay_block.set_dead_time(measurement_delay);



    let initial_temperature = tchx_outlet_temperature_set_point;

    // DRACS hot branch or (mostly) hot leg
    let mut pipe_34 = new_pipe_34(initial_temperature);
    let mut pipe_33 = new_pipe_33(initial_temperature);
    let mut pipe_32 = new_pipe_32(initial_temperature);
    let mut pipe_31a = new_pipe_31a(initial_temperature);
    let mut static_mixer_61_label_31 = new_static_mixer_61_label_31(initial_temperature);
    let mut dhx_tube_side_30b = new_dhx_tube_side_30b(initial_temperature);
    let mut dhx_sthe = new_dhx_sthe_version_1(initial_temperature);
    let mut dhx_tube_side_30a = new_dhx_tube_side_30a(initial_temperature);

    // DRACS cold branch or (mostly) cold leg
    let mut tchx_35a = new_ndhx_tchx_horizontal_35a(initial_temperature);
    let mut tchx_35b = new_ndhx_tchx_vertical_35b(initial_temperature);
    let mut static_mixer_60_label_36 = new_static_mixer_60_label_36(initial_temperature);
    let mut pipe_36a = new_pipe_36a(initial_temperature);
    let mut pipe_37 = new_pipe_37(initial_temperature);
    let mut flowmeter_60_37a = new_flowmeter_60_37a(initial_temperature);
    let mut pipe_38 = new_pipe_38(initial_temperature);
    let mut pipe_39 = new_pipe_39(initial_temperature);

    // pri loop dhx branch top to bottom 5a to 17b 

    let pipe_5a = new_branch_5a(initial_temperature);
    let pipe_26 = new_pipe_26(initial_temperature);
    let pipe_25a = new_pipe_25a(initial_temperature);
    let static_mixer_21_label_25 = new_static_mixer_21_label_25(initial_temperature);
    // here is where the dhx shell side should be (component 24)
    let pipe_23a = new_pipe_23a(initial_temperature);
    let static_mixer_20_label_23 = new_static_mixer_20_label_23(initial_temperature);
    let pipe_22 = new_pipe_22(initial_temperature);
    let flowmeter_20_21a = new_flowmeter_20_label_21a(initial_temperature);
    let pipe_21 = new_pipe_21(initial_temperature);
    let pipe_20 = new_pipe_20(initial_temperature);
    let pipe_19 = new_pipe_19(initial_temperature);
    let pipe_17b = new_branch_17b(initial_temperature);

    // heater branch top to bottom 4 to 18
    let pipe_4 = new_pipe_4(initial_temperature);
    let pipe_3 = new_pipe_3(initial_temperature);
    let pipe_2a = new_pipe_2a(initial_temperature);
    let static_mixer_10_label_2 = new_static_mixer_10_label_2(initial_temperature);
    let heater_top_head_1a = new_heater_top_head_1a(initial_temperature);
    let heater_ver_1 = new_heated_section_version_1_label_1(initial_temperature);
    let heater_bottom_head_1b = new_heater_bottom_head_1b(initial_temperature);
    let pipe_18 = new_pipe_18(initial_temperature);




    todo!()

}
