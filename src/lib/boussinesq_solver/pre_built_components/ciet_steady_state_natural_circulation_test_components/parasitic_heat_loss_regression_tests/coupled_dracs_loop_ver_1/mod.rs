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
pub fn regression_long_test_uncalibrated_dracs_loop(){

    use std::thread;

    use regression_coupled_dracs_loop_version_1::*;



    let thread_1 = thread::spawn(||{
        let max_simulation_time_seconds: f64 = 3000.0;
        // expect overprediction of mass flowrates in both loops 
        // to about 10%
        let pri_loop_relative_tolerance = 0.1;
        let dracs_loop_relative_tolerance = 0.1;

        // I'm writing in this format so that the data will be easier 
        // to copy over to csv
        let (heater_power_watts,
            tchx_outlet_temp_degc,
            experimental_primary_mass_flowrate_kg_per_s,
            experimental_dracs_mass_flowrate_kg_per_s,
            simulated_expected_primary_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s) 
            = (2764.53, 40.0, 4.6990e-2, 3.5470e-2, 4.9890e-2, 3.8766e-2);

        dbg!(max_simulation_time_seconds,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance);


        regression_coupled_dracs_loop_version_1(
            heater_power_watts, 
            max_simulation_time_seconds,
            tchx_outlet_temp_degc,
            experimental_primary_mass_flowrate_kg_per_s,
            experimental_dracs_mass_flowrate_kg_per_s,
            simulated_expected_primary_mass_flowrate_kg_per_s,
            simulated_expected_dracs_mass_flowrate_kg_per_s,
            pri_loop_relative_tolerance,
            dracs_loop_relative_tolerance,
        ).unwrap();

    });


    thread_1.join().unwrap();
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

