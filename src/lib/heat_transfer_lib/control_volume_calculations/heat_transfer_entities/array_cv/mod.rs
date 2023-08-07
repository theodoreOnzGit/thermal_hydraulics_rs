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


use std::ops::{Add, Mul, Sub, Div, Rem};

use ndarray::*;
use ndarray_linalg::*;
use uom::num_traits::Num;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermal_conductance::watt_per_kelvin;


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

    let mut M: Array2<ThermalConductance> = Array::zeros((3,3));

    M[[0,0]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    M[[0,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    M[[0,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    M[[1,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    M[[1,1]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    M[[1,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    M[[2,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    M[[2,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    M[[2,2]] = ThermalConductance::new::<watt_per_kelvin>(3.0);

    // this is the power vector
    let mut S: Array1<Power> = Array::zeros(3);

    S[0] = Power::new::<watt>(1.0);
    S[1] = Power::new::<watt>(2.0);
    S[2] = Power::new::<watt>(3.01);


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

    Ok(())
}


#[test] 
fn test_linalg(){
    solve_example_array().unwrap();
    factorize().unwrap();

    sandbox().unwrap();
}

///// I created an internal version of the power struct 
///// that specifically does okay for linalg calculations
//#[derive(Clone,Copy,Debug, PartialOrd, PartialEq)]
//pub struct LinalgPower {
//    power: Power
//}
//
//impl cauchy::Scalar for LinalgPower {}
//impl Num for LinalgPower {
//    type FromStrRadixErr;
//
//    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
//        todo!()
//    }
//}
//
//
////impl LinalgScalar for LinalgPower {
////    
////}
////
//impl uom::num_traits::One for LinalgPower {
//    fn one() -> Self {
//        // gives a power of one watt 
//        //
//
//        let power = Power::new::<watt>(1.0);
//        return LinalgPower { power };
//    }
//}
//
//impl uom::num_traits::Zero for LinalgPower {
//    fn zero() -> Self {
//        let power = Power::new::<watt>(0.0);
//        return LinalgPower { power };
//    }
//
//    fn is_zero(&self) -> bool {
//        let power_to_test: Power = self.power;
//
//        // convert to f64 and return the answer
//        return power_to_test.value.is_zero();
//    }
//}
//
//impl Add for LinalgPower{}
//
//impl Mul for LinalgPower{}
//
//impl Sub for LinalgPower{}
//
//impl Div for LinalgPower{
//    type Output = LinalgPower;
//
//    fn div(self, rhs: Self) -> Self::Output {
//
//        let power: Power = self.power;
//        let rhs: Power = rhs.power;
//
//        let div: Power = power.div(rhs);
//        return Self { power:div };
//
//        // problem, if i do it for power, I need to do it for all 
//        // quantities =(
//    }
//}
//
//impl Rem for LinalgPower {
//    type Output = LinalgPower;
//
//    fn rem(self, rhs: Self) -> Self::Output {
//        // first get the power 
//
//        let power: Power = self.power;
//        let rhs: Power = rhs.power;
//
//        let rem: Power = power.rem(rhs);
//        return Self { power:rem };
//    }
//}
