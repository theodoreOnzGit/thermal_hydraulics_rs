use uom::si::{f64::*, mass_rate::kilogram_per_second, power::watt};

use crate::prelude::beta_testing::ThermalHydraulicsLibError;
#[test]
pub fn case_a(){


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


#[test]
pub fn case_b(){

    verify_isolated_dhx_analytical_solution(
        454.4, 
        2.4297e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        766.2, 
        3.0478e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1004.4, 
        3.4263e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1211.2, 
        3.6997e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1409.0, 
        3.9461e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1607.4, 
        4.1702e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1804.6, 
        4.3738e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        2004.9, 
        4.5643e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        2211.0, 
        4.7497e-2
        ).unwrap();
}

#[test]
pub fn case_c(){

    verify_isolated_dhx_analytical_solution(
        582.6, 
        2.7989e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        785.9, 
        3.1748e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        971.4, 
        3.4616e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1185.2, 
        3.7682e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1369.1, 
        4.0000e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1584.1, 
        4.2382e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1763.7, 
        4.4318e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        1970.0, 
        4.6319e-2
        ).unwrap();
    verify_isolated_dhx_analytical_solution(
        2177.0, 
        4.8155e-2
        ).unwrap();
}

/// function to verify the dhx analytical solution
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

/// debugging tests for functions to make natural circulation 
/// testing easier 
pub mod debugging;
