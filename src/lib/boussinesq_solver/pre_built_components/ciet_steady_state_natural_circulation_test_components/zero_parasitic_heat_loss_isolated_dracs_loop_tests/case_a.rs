
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
pub fn case_a_tchx_out_319_kelvin_46_celsius(){
    use uom::si::{f64::*, mass_rate::kilogram_per_second, power::watt};

    use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

    fn verify_isolated_dhx_analytical_solution(
        input_power_watts: f64,
        analytical_solution_mass_flowrate_kg_per_s: f64) -> 
        Result<(),ThermalHydraulicsLibError>{

            let input_power = Power::new::<watt>(input_power_watts);
            let analytical_solution_mass_flowrate = 
                MassRate::new::<kilogram_per_second>(
                    analytical_solution_mass_flowrate_kg_per_s);

            // max error is 0.5% according to SAM 
            // is okay, because typical flowmeter measurement error is 2% anyway
            todo!()

    }

    verify_isolated_dhx_analytical_solution(
        931.8, 
        3.4967e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1088.3, 
        3.7214e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1338.4, 
        4.0525e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1470.6, 
        4.2045e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1699.9, 
        4.4583e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1876.5, 
        4.6309e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        2137.0, 
        4.8754e-2
        ).unwrap();
}
