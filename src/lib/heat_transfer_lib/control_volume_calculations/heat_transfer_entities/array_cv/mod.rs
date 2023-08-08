// to start off array_cv, we take inspiration from GeN-Foam's 
// lumpedNuclearStructure code 
//
// lumpedNuclearStructure was meant to model a heat generating 
// pebble geometry because GeN-Foam used to only allow 
// for pin shaped geometry. To deal with this, we have several 
// control volumes or nodes with thermal conductances between 
// the nodes
//
//
//
// now, for matrix solution, i use intel-mkl-static in ndarray_linalg 
// library 
//
// this is because ndarray_linalg using intel-mkl-static is cross 
// platform, and it can be used for windows, macos and linux 
//
// secondly, the intel-mkl-static library compiles the lapack library 
// locally and links it statically, rather than at the system level 
// therefore, the user won't have to worry as much about system 
// dependencies which can cause some headache
//
// btw, in future implementations of thermal_hydraulics_rs, i might 
// want to copy and paste how ndarray-linalg constructs its error types 
// so that my results are similar



use ndarray::*;
use ndarray_linalg::*;
use uom::ConstZero;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::thermal_conductance::watt_per_kelvin;

/// This module contains direct translations of GeN-Foam code 
/// into rust for the lumped nuclear structure
#[path = "./gen-foam-lumped-nuclear-structure.rs"]
pub mod gen_foam_lumped_nuclear_structure;
use gen_foam_lumped_nuclear_structure::solve_conductance_matrix_power_vector;
/// This module contains adapts the GeN-Foam code for the 
/// lumped nuclear structure
/// for the heat transfer module
pub mod lumped_nuclear_structure_inspired_functions;
use gen_foam_lumped_nuclear_structure::*;

/// this module contains constructors for an array cv which 
/// functions as a one dimensional cartesian conduction medium
pub mod one_dimension_cartesian_conducting_medium;
pub use one_dimension_cartesian_conducting_medium::*;


// Solve `Ax=b`
fn solve_example_array() -> Result<(), error::LinalgError> {
    let a: Array2<f64> = 
    arr2(&[[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]]);
    let b: Array1<f64> = arr1(&[1.0,2.0,3.01]);
    let _x = a.solve(&b)?;
    Ok(())
}

// Solve `Ax=b` for many b with fixed A
fn factorize() -> Result<(), error::LinalgError> {
    let a: Array2<f64> = random((3, 3));
    let f = a.factorize_into()?; // LU factorize A (A is consumed)
    for _ in 0..10 {
        let b: Array1<f64> = random(3);
        let _x = f.solve_into(b)?; // solve Ax=b using factorized L, U
    }
    Ok(())
}

// sandbox 
fn sandbox() -> Result<(), error::LinalgError> {

    // i want to index into the array and change stuff
    let mut a: Array2<f64> = 
    arr2(&[[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]]);


    // https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html
    //
    // the syntax to access by element is slightly different from c
    a[[1,2]] = 0.5;



    let b: Array1<f64> = arr1(&[1.0,2.0,3.01]);
    let _x = a.solve(&b)?;


    // we can also initialise an array using 
    // 
    // surprisingly, we can use thermal conductance as well
    // and hopefully in general, 

    let mut thermal_conductance_matrix: Array2<ThermalConductance> = Array::zeros((3,3));

    thermal_conductance_matrix[[0,0]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    thermal_conductance_matrix[[0,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[0,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    thermal_conductance_matrix[[1,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[1,1]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    thermal_conductance_matrix[[1,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    thermal_conductance_matrix[[2,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[2,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[2,2]] = ThermalConductance::new::<watt_per_kelvin>(3.0);

    // this is the power vector
    let mut power_vector: Array1<Power> = Array::zeros(3);

    power_vector[0] = Power::new::<watt>(1.0);
    power_vector[1] = Power::new::<watt>(2.0);
    power_vector[2] = Power::new::<watt>(3.01);


    // solve linear system
    let _temperature_vector: Array1<ThermodynamicTemperature> = 
    solve_conductance_matrix_power_vector(thermal_conductance_matrix,
        power_vector)?;


    // now how to solve this? I'll probably have to convert this into 
    // a float, but I lose my unit safety
    //
    // we can do type conversion here: 
    // https://github.com/rust-ndarray/ndarray/blob/master/examples/type_conversion.rs
    //
    // direct solving is problematic

    //let T = M.solve(&S)?;
    //
    // this is because ndarray allows you to create non dimensional 
    // arrays of anytime but ndarray-linalg only deals with f64
    // or other numeric types
    //
    // one other way is to convert all into f64, and verify the solution 
    // by mapping it back to its appropriate units 
    //
    // so T will be an f64 vector, which should be temperature 
    // S will be an f64 vector as well 
    // M will be an f64 matrix 
    //
    // so in effect, a custom solve method with appropriate unit 
    // checks would suffice
    //

    Ok(())
}
