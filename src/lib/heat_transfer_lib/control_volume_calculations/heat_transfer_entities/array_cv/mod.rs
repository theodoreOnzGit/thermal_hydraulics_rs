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
use uom::ConstZero;
use uom::num_traits::Num;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;
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


#[inline]
fn solve_conductance_matrix_power_vector(
    thermal_conductance_matrix: Array2<ThermalConductance>,
    power_vector: Array1<Power>)
-> Result<Array1<ThermodynamicTemperature>, error::LinalgError>{

    // I can of course convert it into f64 types 
    //
    //

    let get_value_conductance = |conductance: &ThermalConductance| {
        return conductance.value;
    };

    let get_value_power = |power: &Power| {
        return power.value;
    };

    // i'm allowing non snake case so that the syntax is the same as 
    // GeN-Foam
    #[allow(non_snake_case)]
    let M: Array2<f64> = 
    thermal_conductance_matrix.map(get_value_conductance);

    #[allow(non_snake_case)]
    let S: Array1<f64> = power_vector.map(get_value_power);

    // now for the raw temperature matrix 

    #[allow(non_snake_case)]
    let T: Array1<f64> = M.solve(&S)?;

    // To check for unit safety, I can just perform one calc

    let _unit_check: Power = 
    power_vector[0] + 
    thermal_conductance_matrix[[0,0]] 
    * ThermodynamicTemperature::ZERO;

    // now map T back to a ThermodynamicTemperature
    // T is already a ThermodynamicTemperature, so don't need to manually 
    // convert, do it in kelvin 
    //

    let convert_f64_to_kelvin_temperature = |float: &f64| {
        return ThermodynamicTemperature::new::<kelvin>(*float);
    };


    // this is the last step
    let temperature_vector: Array1<ThermodynamicTemperature> 
    = T.map(convert_f64_to_kelvin_temperature);

    return Ok(temperature_vector);
}


#[test] 
fn test_linalg(){
    solve_example_array().unwrap();
    factorize().unwrap();

    sandbox().unwrap();
}

// can't do this
//impl<D, U, V> LinalgScalar for uom::si::Quantity<D, U, V> {}

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
