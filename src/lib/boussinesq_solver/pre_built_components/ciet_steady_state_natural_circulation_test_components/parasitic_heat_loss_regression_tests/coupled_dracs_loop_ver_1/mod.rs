#[cfg(test)]
#[test]
pub fn quick_test_uncalibrated_dracs_loop(){

    use validate_coupled_dracs_loop_version_1::*;
    let max_simulation_time_seconds: f64 = 400.0;
    let pri_loop_relative_tolerance = 0.03;
    let dracs_loop_relative_tolerance = 0.4;

    validate_coupled_dracs_loop_version_1(
        2764.53, 
        max_simulation_time_seconds,
        40.0,
        4.6990e-2,
        3.5470e-2,
        pri_loop_relative_tolerance,
        dracs_loop_relative_tolerance,
        ).unwrap();
}
#[cfg(test)]
#[test]
pub fn long_test_uncalibrated_dracs_loop(){

    use validate_coupled_dracs_loop_version_1::*;
    let max_simulation_time_seconds: f64 = 4000.0;
    // expect overprediction of mass flowrates in both loops 
    // to about 10%
    let pri_loop_relative_tolerance = 0.1;
    let dracs_loop_relative_tolerance = 0.1;

    validate_coupled_dracs_loop_version_1(
        2764.53, 
        max_simulation_time_seconds,
        40.0,
        4.6990e-2,
        3.5470e-2,
        pri_loop_relative_tolerance,
        dracs_loop_relative_tolerance,
        ).unwrap();
}
#[cfg(test)]
#[test]
pub fn regression_long_test_uncalibrated_dracs_loop_set_c(){

    use std::thread;

    use regression_coupled_dracs_loop_version_1::*;

    let set_c1 = thread::Builder::new()
        .name("set_c1".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (841.02, 40.0, 2.6860e-2, 2.0030e-2, 3.0218e-2, 2.2822e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c2 = thread::
        Builder::new().name("set_c2".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (1158.69, 40.0, 3.0550e-2, 2.3670e-2, 3.4997e-2, 2.6789e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();
    let set_c3 = thread::
        Builder::new().name("set_c3".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (1409.22, 40.0, 3.3450e-2, 2.6350e-2, 3.8177e-2, 2.9411e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c4 = thread::
        Builder::new().name("set_c4".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (1736.11, 40.0, 3.6490e-2, 2.9490e-2, 4.1805e-2, 3.2234e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c5 = thread::
        Builder::new().name("set_c5".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (2026.29, 40.0, 3.8690e-2, 3.1900e-2, 4.4660e-2, 3.4584e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c6 = thread::
        Builder::new().name("set_c6".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (2288.83, 40.0, 4.1150e-2, 3.4120e-2, 4.7011e-2, 3.6355e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c7 = thread::
        Builder::new().name("set_c7".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (2508.71, 40.0, 4.3120e-2, 3.5620e-2, 4.8844e-2, 3.7680e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c8 = thread::
        Builder::new().name("set_c8".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 15%
        let pri_loop_relative_tolerance = 0.15;
        let dracs_loop_relative_tolerance = 0.15;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (2685.83, 40.0, 4.5090e-2, 3.5930e-2, 5.0242e-2, 3.8670e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();

    let set_c9 = thread::
        Builder::new().name("set_c9".to_string()).spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 10%
        let pri_loop_relative_tolerance = 0.1;
        let dracs_loop_relative_tolerance = 0.1;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s) 
            = (2764.53, 40.0, 4.6990e-2, 3.5470e-2, 5.0843e-2, 3.9091e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_dracs_mass_flowrate_kg_per_s,
            experimental_pri_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            simulated_expected_pri_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    }).unwrap();


    set_c1.join().unwrap();
    set_c2.join().unwrap();
    set_c3.join().unwrap();
    set_c4.join().unwrap();
    set_c5.join().unwrap();
    set_c6.join().unwrap();
    set_c7.join().unwrap();
    set_c8.join().unwrap();
    set_c9.join().unwrap();
}


/// function to test uncalibrated 
/// coupled dracs loop and compare with experimental data 
/// this is more of a regression function, so I want to check the 
/// output of the uncalibrated loop
///
/// the DHX here uses uncalibrated Gnielinski correlations 
/// to estimate heat transfer coefficients
pub mod regression_coupled_dracs_loop_version_1;


/// function to validate coupled DRACS loop to experimental data 
/// within a given tolerance
/// version 1,
/// the DHX here uses uncalibrated Gnielinski correlations 
/// to estimate heat transfer coefficients
pub mod validate_coupled_dracs_loop_version_1;

