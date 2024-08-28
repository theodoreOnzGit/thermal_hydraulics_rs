use uom::si::{f64::*, mass_rate::kilogram_per_second, power::watt};

use crate::prelude::beta_testing::ThermalHydraulicsLibError;

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

    verify_coupled_dhx_analytical_solution(
        1479.86, 
        3.3410e-2,
        2.7380e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1653.90, 
        3.5440e-2,
        2.8190e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2014.51, 
        3.8770e-2,
        3.2360e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2178.49, 
        4.0110e-2,
        3.2550e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2395.90, 
        4.2770e-2,
        3.3900e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2491.87, 
        4.4650e-2,
        3.3550e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2696.24, 
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

    verify_coupled_dhx_analytical_solution(
        655.16, 
        2.3290e-2,
        1.7310e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1054.32, 
        2.9520e-2,
        2.1980e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1394.70, 
        3.3240e-2,
        2.5700e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1685.62, 
        3.6110e-2,
        2.8460e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1987.75, 
        3.8410e-2,
        3.1180e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2282.01, 
        4.0630e-2,
        3.3740e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2546.60, 
        4.2700e-2,
        3.5770e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2874.03, 
        4.4560e-2,
        3.7960e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        3031.16, 
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

    verify_coupled_dhx_analytical_solution(
        841.02, 
        2.6860e-2,
        2.0030e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1158.69, 
        3.0550e-2,
        2.3670e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1409.22, 
        3.3450e-2,
        2.6350e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        1736.11, 
        3.6490e-2,
        2.9490e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2026.29, 
        3.8690e-2,
        3.1900e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2288.83, 
        4.1150e-2,
        3.4120e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2508.71, 
        4.3120e-2,
        3.5620e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2685.83, 
        4.5090e-2,
        3.5930e-2,
        ).unwrap();
    verify_coupled_dhx_analytical_solution(
        2764.53, 
        4.6990e-2,
        3.4570e-2,
        ).unwrap();
}

/// function to verify the dhx analytical solution
pub fn verify_coupled_dhx_analytical_solution(
    input_power_watts: f64,
    experimental_dracs_mass_flowrate_kg_per_s: f64,
    experimental_primary_mass_flowrate_kg_per_s: f64) -> 
Result<(),ThermalHydraulicsLibError>{

    let _input_power = Power::new::<watt>(input_power_watts);
    let _experimental_dracs_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_dracs_mass_flowrate_kg_per_s);
    let _experimental_primary_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_primary_mass_flowrate_kg_per_s);

    // max error is 0.5% according to SAM 
    // is okay, because typical flowmeter measurement error is 2% anyway
    todo!()

}

/// debugging tests for functions to make natural circulation 
/// testing easier 
pub mod debugging;
