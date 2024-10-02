#![warn(missing_docs)]
extern crate peroxide;
use peroxide::prelude::*;

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

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
pub fn custom_f_ldk(custom_darcy: &dyn Fn(f64, f64) -> Result<f64,ThermalHydraulicsLibError>,
        reynolds_number: f64,
        roughness_ratio: f64,
        length_to_diameter_ratio: f64,
        custom_k: &dyn Fn(f64) -> Result<f64,ThermalHydraulicsLibError>) 
    -> Result<f64,ThermalHydraulicsLibError>{

    if roughness_ratio < 0.0 {
        panic!("roughnessRatio<0.0");
    }

    if length_to_diameter_ratio <= 0.0 {
        panic!("lengthToDiameterRatio<=0.0");
    }

    let k = custom_k(reynolds_number)?;


    let f = custom_darcy(reynolds_number, roughness_ratio)?;
    let f_ldk = f*length_to_diameter_ratio + k;

    return Ok(f_ldk);
}

/// this is a special case of the fLDK component,
/// where we just specify a custom K but friction factor is based
/// on darcy friction factor
pub fn custom_kpipe(reynolds_number: f64,
    roughness_ratio: f64,
    length_to_diameter_ratio: f64,
    custom_k: &dyn Fn(f64) -> Result<f64,ThermalHydraulicsLibError>) -> 
Result<f64, ThermalHydraulicsLibError>{

    let darcy_fn = crate::boussinesq_solver::fluid_mechanics_correlations::
        darcy;

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
pub fn custom_kpipe_be_d(reynolds_number: f64,
                    roughness_ratio: f64,
                    length_to_diameter_ratio: f64,
                    custom_k: &dyn Fn(f64) -> Result<f64,ThermalHydraulicsLibError>) 
    -> Result<f64,ThermalHydraulicsLibError>{

    if reynolds_number == 0.0 {
        return Ok(0.0);
    }

    let f_ldk = custom_kpipe(reynolds_number,
                           roughness_ratio,
                           length_to_diameter_ratio,
                           custom_k)?;

    let bejan_d = 0.5*f_ldk*reynolds_number.powf(2.0);

    return Ok(bejan_d);

}


/// this functions calculates the bejan number using the
/// custom fLDK formula
pub fn custom_f_ldk_be_d(custom_darcy: &dyn Fn(f64, f64) -> Result<f64,ThermalHydraulicsLibError>, 
                        reynolds_number: f64,
                        roughness_ratio: f64,
                        length_to_diameter_ratio: f64,
                        custom_k: &dyn Fn(f64) -> 
                        Result<f64,ThermalHydraulicsLibError>) -> Result<f64,ThermalHydraulicsLibError>{

    if reynolds_number == 0.0 {
        return Ok(0.0);
    }

    let f_ldk = custom_f_ldk(custom_darcy,
                           reynolds_number,
                           roughness_ratio,
                           length_to_diameter_ratio,
                           custom_k)?;

    let bejan_d = 0.5*f_ldk*reynolds_number.powf(2.0);

    return Ok(bejan_d);

}

/// this code allos us to get Reynold's number from a Bejan
/// number for a custom pipe.
/// i make no assumptions about the symmetry of flow
/// ie. i don't make assumptions about whether
/// the pipe exhibits the same pressure loss
/// in forwards and backwards flow,
/// that is up to the user to decide when 
/// customDarcy and customK is put in
pub fn get_reynolds(
    custom_darcy: &dyn Fn(f64, f64) -> Result<f64,ThermalHydraulicsLibError>, 
    bejan_d: f64,
    roughness_ratio: f64,
    length_to_diameter: f64,
    custom_k: &dyn Fn(f64) -> Result<f64,ThermalHydraulicsLibError>) -> Result<f64,ThermalHydraulicsLibError> {

    if length_to_diameter <= 0.0 {
        panic!("lengthToDiameterRatio<=0.0");
    }

    if roughness_ratio < 0.0 {
        panic!("roughnessRatio<0.0");
    }


    // this part deals with negative Be_L values
    // invalid Be_L values

    let max_reynolds = 1.0e12;

    // i calculate the Be_D corresponding to 
    // Re = 1e12
    let max_bejan_d = custom_f_ldk_be_d(
        custom_darcy,
        max_reynolds,
        roughness_ratio, 
        length_to_diameter,
        custom_k)?;

    if bejan_d >= max_bejan_d {
        panic!("Be too large");
    }
    // the above checks for all the relevant exceptions
    // including formLossK < 0
    //
    // now we are ready to do root finding
    //
    // the underlying equation is 
    // Be = 0.5*fLDK*Re^2


    let pressure_drop_root = |reynolds: AD| -> AD {
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

        let reynolds_double = reynolds.x();
        let f_ldk_term = custom_f_ldk_be_d(
            custom_darcy,
            reynolds_double, 
            roughness_ratio,
            length_to_diameter,
            custom_k).unwrap();

        return AD0(bejan_d - f_ldk_term);

    };

    let reynolds_number_result = bisection(pressure_drop_root,
        (-max_reynolds,max_reynolds),
        100,
        1e-8);



    // the unwrap turns the result into f64
    let reynolds_number = reynolds_number_result.unwrap();

    return Ok(reynolds_number);
}


