
/// functions used for calculating the thermal hydraulics inside the DRACS 
/// loop
pub mod dracs_loop_calc_functions;

/// functions used for calculating the thermal hydraulics inside 
/// the Heater and DHX branch 
/// Note: heater v1.0 is used
pub mod pri_loop_calc_functions;

#[test]
/// We use:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 4 provides the data we use here
/// 
///
#[ignore = "under construction"]
pub fn case_a_tchx_out_319_kelvin_46_celsius(){

    // data is:
    //
    // 1) heat added in watts
    // 2) experimental DRACS flowrate in kg/s
    // 3) experimental Primary flowrate in kg/s

    verify_coupled_dhx_analytical_solution_version_1(
        1479.86, 
        46.0,
        3.3410e-2,
        2.7380e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1653.90, 
        46.0,
        3.5440e-2,
        2.8190e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2014.51, 
        46.0,
        3.8770e-2,
        3.2360e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2178.49, 
        46.0,
        4.0110e-2,
        3.2550e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2395.90, 
        46.0,
        4.2770e-2,
        3.3900e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2491.87, 
        46.0,
        4.4650e-2,
        3.3550e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2696.24, 
        46.0,
        4.7100e-2,
        3.4620e-2,
        ).unwrap();
}


/// In the original SAM publication
///
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// I found it hard to distinguish what TCHX temperatures case A,
/// B and C were.
///
/// But there was another publication which shows which is test group 
/// corresponds to which temperature:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 3 also provides the data 
/// 
///
#[test]
#[ignore = "under construction"]
pub fn case_b_tchx_out_308_kelvin_35_celsius(){

    verify_coupled_dhx_analytical_solution_version_1(
        655.16, 
        35.0,
        2.3290e-2,
        1.7310e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1054.32, 
        35.0,
        2.9520e-2,
        2.1980e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1394.70, 
        35.0,
        3.3240e-2,
        2.5700e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1685.62, 
        35.0,
        3.6110e-2,
        2.8460e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1987.75, 
        35.0,
        3.8410e-2,
        3.1180e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2282.01, 
        35.0,
        4.0630e-2,
        3.3740e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2546.60, 
        35.0,
        4.2700e-2,
        3.5770e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2874.03, 
        35.0,
        4.4560e-2,
        3.7960e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        3031.16, 
        35.0,
        4.6360e-2,
        3.8490e-2,
        ).unwrap();
}

/// In the original SAM publication
///
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// I found it hard to distinguish what TCHX temperatures case A,
/// B and C were.
///
/// But there was another publication which shows which is test group 
/// corresponds to which temperature:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 3 also provides the data 
/// 
///
#[test]
#[ignore = "under construction"]
pub fn case_c_tchx_out_313_kelvin_40_celsius(){

    verify_coupled_dhx_analytical_solution_version_1(
        841.02, 
        40.0,
        2.6860e-2,
        2.0030e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1158.69, 
        40.0,
        3.0550e-2,
        2.3670e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1409.22, 
        40.0,
        3.3450e-2,
        2.6350e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        1736.11, 
        40.0,
        3.6490e-2,
        2.9490e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2026.29, 
        40.0,
        3.8690e-2,
        3.1900e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2288.83, 
        40.0,
        4.1150e-2,
        3.4120e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2508.71, 
        40.0,
        4.3120e-2,
        3.5620e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2685.83, 
        40.0,
        4.5090e-2,
        3.5930e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution_version_1(
        2764.53, 
        40.0,
        4.6990e-2,
        3.4570e-2,
        ).unwrap();
}

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
    let mut dhx_tube_side_heat_exchanger_30 = new_isolated_dhx_tube_side_30(initial_temperature);
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
    todo!()

}
/// constructor for the dhx shell and tube heat exchanger 
/// based on Zou's specifications
pub mod dhx_constructor;

/// debugging tests for functions to make natural circulation 
/// testing easier 
pub mod debugging;
