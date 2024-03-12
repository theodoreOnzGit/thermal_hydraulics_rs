#![warn(missing_docs)]
extern crate peroxide;
use peroxide::prelude::*;

use crate::prelude::alpha_nightly::ThermalHydraulicsLibError;

// This library was developed for use in my PhD thesis under supervision 
// of Professor Per F. Peterson. It is part of a thermal hydraulics
// library in Rust that is released under the GNU General Public License
// v 3.0. This is partly due to the fact that some of the libraries 
// inherit from GeN-Foam and OpenFOAM, both licensed under GNU General
// Public License v3.0.
//
// As such, the entire library is released under GNU GPL v3.0. It is a strong 
// copyleft license which means you cannot use it in proprietary software.
//
//
// License
//    This is file is part of a thermal hydraulics library written 
//    in rust meant to help with the
//    fluid mechanics and heat transfer aspects of the calculations
//    for the Compact Integral Effects Tests (CIET) and hopefully 
//    Gen IV Reactors such as the Fluoride Salt cooled High Temperature 
//    Reactor (FHR)
//     
//    Copyright (C) 2022-2023  Theodore Kay Chen Ong, Singapore Nuclear
//    Research and Safety Initiative, Per F. Peterson, University of 
//    California, Berkeley Thermal Hydraulics Laboratory
//
//    thermal_hydrualics_rs is free software; you can 
//    redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by the
//    Free Software Foundation; either version 2 of the License, or (at your
//    option) any later version.
//
//    thermal_hydrualics_rs is distributed in the hope 
//    that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    This thermal hydraulics library 
//    contains some code copied from GeN-Foam, and OpenFOAM derivative.
//    This offering is not approved or endorsed by the OpenFOAM Foundation nor
//    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
//    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.
//    Nor is it endorsed by the authors and owners of GeN-Foam.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// © All rights reserved. Theodore Kay Chen Ong,
// Singapore Nuclear Research and Safety Initiative,
// Per F. Peterson,
// University of California, Berkeley Thermal Hydraulics Laboratory
//
// Main author of the code: Theodore Kay Chen Ong, supervised by
// Professor Per F. Peterson
//
// Btw, I have no affiliation with the Rust Foundation.

// here are the functions used for components with 
// custom friction factor and K, rather messy but


/// this first function allows for custom fldk, 
/// ie both friction factor and form loss k are user defined
/// <https://stackoverflow.com/questions/36390665/how-do-you-pass-a-rust-function-as-a-parameter>
pub fn custom_f_ldk(custom_darcy: &dyn Fn(f64, f64) -> f64,
        reynolds_number: f64,
        roughness_ratio: f64,
        length_to_diameter_ratio: f64,
        custom_k: &dyn Fn(f64) -> f64) -> Result<f64,ThermalHydraulicsLibError>{

    if roughness_ratio < 0.0 {
        panic!("roughnessRatio<0.0");
    }

    if length_to_diameter_ratio <= 0.0 {
        panic!("lengthToDiameterRatio<=0.0");
    }

    let k = custom_k(reynolds_number);


    let f = custom_darcy(reynolds_number, roughness_ratio);
    let f_ldk = f*length_to_diameter_ratio + k;

    return Ok(f_ldk);
}

/// this is a special case of the fLDK component,
/// where we just specify a custom K but friction factor is based
/// on darcy friction factor
pub fn custom_Kpipe(reynolds_number: f64,
    roughness_ratio: f64,
    length_to_diameter_ratio: f64,
    custom_k: &dyn Fn(f64) -> f64) -> 
Result<f64, ThermalHydraulicsLibError>{

    let darcy_fn = crate::fluid_mechanics_lib::
        churchill_friction_factor::darcy;

    let f_ldk_result = custom_f_ldk(&darcy_fn,
        reynolds_number,
        roughness_ratio,
        length_to_diameter_ratio,
        custom_k);

    return f_ldk_result;

}


/// This function is special,
/// not really used, as often
/// it assumes that the form losses, K for the pipe 
/// take some functional form rather than staying constant
///
#[allow(non_snake_case)]
pub fn custom_Kpipe_Be_D(ReynoldsNumber: f64,
                    roughnessRatio: f64,
                    lengthToDiameterRatio: f64,
                    customK: &dyn Fn(f64) -> f64) -> Result<f64,ThermalHydraulicsLibError>{

    if ReynoldsNumber == 0.0 {
        return Ok(0.0);
    }

    let fLDK = custom_Kpipe(ReynoldsNumber,
                           roughnessRatio,
                           lengthToDiameterRatio,
                           customK)?;

    let Be_D = 0.5*fLDK*ReynoldsNumber.powf(2.0);

    return Ok(Be_D);

}


#[allow(non_snake_case)]
/// this functions calculates the bejan number using the
/// custom fLDK formula
pub fn custom_fLDK_Be_D(customDarcy: &dyn Fn(f64, f64) -> f64, 
                        ReynoldsNumber: f64,
                        roughnessRatio: f64,
                        lengthToDiameterRatio: f64,
                        customK: &dyn Fn(f64) -> f64) -> Result<f64,ThermalHydraulicsLibError>{

    if ReynoldsNumber == 0.0 {
        return Ok(0.0);
    }

    let fLDK = custom_f_ldk(customDarcy,
                           ReynoldsNumber,
                           roughnessRatio,
                           lengthToDiameterRatio,
                           customK)?;

    let Be_D = 0.5*fLDK*ReynoldsNumber.powf(2.0);

    return Ok(Be_D);

}

/// this code allos us to get Reynold's number from a Bejan
/// number for a custom pipe.
/// i make no assumptions about the symmetry of flow
/// ie. i don't make assumptions about whether
/// the pipe exhibits the same pressure loss
/// in forwards and backwards flow,
/// that is up to the user to decide when 
/// customDarcy and customK is put in
#[allow(non_snake_case)]
pub fn getRe(customDarcy: &dyn Fn(f64, f64) -> f64, 
             Be_D: f64,
             roughnessRatio: f64,
             lengthToDiameter: f64,
             customK: &dyn Fn(f64) -> f64) -> Result<f64,ThermalHydraulicsLibError> {

    if lengthToDiameter <= 0.0 {
        panic!("lengthToDiameterRatio<=0.0");
    }

    if roughnessRatio < 0.0 {
        panic!("roughnessRatio<0.0");
    }


    // this part deals with negative Be_L values
    // invalid Be_L values

    let maxRe = 1.0e12;

    // i calculate the Be_D corresponding to 
    // Re = 1e12
    let maxBe_D = custom_fLDK_Be_D(
        customDarcy,
        maxRe,
        roughnessRatio, 
        lengthToDiameter,
        customK)?;

    if Be_D >= maxBe_D {
        panic!("Be too large");
    }
    // the above checks for all the relevant exceptions
    // including formLossK < 0
    //
    // now we are ready to do root finding
    //
    // the underlying equation is 
    // Be = 0.5*fLDK*Re^2


    let pressureDropRoot = |Re: AD| -> AD {
        // i'm solving for
        // Be - 0.5*fLDK*Re^2 = 0 
        // the fLDK term can be calculated using
        // getBe
        //
        // now i don't really need the interpolation
        // term in here because when Re = 0,
        // Be = 0 in the getBe code.
        // so really, no need for fancy interpolation.
        //
        // Now in peroxide, the type taken in and out
        // is not a f64 double
        // but rather AD which stands for automatic 
        // differentiation
        // https://docs.rs/peroxide/latest/peroxide/structure/ad/index.html

        let reynoldsDouble = Re.x();
        let fLDKterm = custom_fLDK_Be_D(
            customDarcy,
            reynoldsDouble, 
            roughnessRatio,
            lengthToDiameter,
            customK).unwrap();

        return AD0(Be_D - fLDKterm);

    };

    let ReynoldsNumberResult = bisection(pressureDropRoot,
                                         (-maxRe,maxRe),
                                         100,
                                         1e-8);



    // the unwrap turns the result into f64
    let ReynoldsNumber = ReynoldsNumberResult.unwrap();

    return Ok(ReynoldsNumber);
}


