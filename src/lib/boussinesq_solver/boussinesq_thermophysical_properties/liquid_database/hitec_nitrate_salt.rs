#[warn(missing_docs)]

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
//    Copyright (C) 2022-2024  Theodore Kay Chen Ong, Singapore Nuclear
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
// Btw, I have no affiliation with the Rust foundation.
use uom::si::f64::*;
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::dynamic_viscosity::{millipascal_second, pascal_second};
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::available_energy::joule_per_kilogram;

// this is for the root finding algorithms
extern crate peroxide;
use peroxide::prelude::*;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::{range_check, LiquidMaterial, Material};
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

/// function to obtain nitrate salt density
/// given a temperature
///
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68./// Jana, S. S., 
/// Maheshwari, N. K., & Vijayan, P. K. (2016). 
/// 
///
/// rho (kg/m3) = 2280.22  - 0.773 T(K)
pub fn get_hitec_equimolar_density(
    fluid_temp: ThermodynamicTemperature) -> Result<MassDensity,ThermalHydraulicsLibError> {


    // first we check if fluid temp is between 440-800 K (range of validity)
    // panic otherwise
    range_check_hitec_salt(fluid_temp)?;
    let fluid_temp_kelvin = fluid_temp.get::<kelvin>();
    let a = 2280.22;
    let b = -0.773;
    // generic correlation is:
    // a + bT + cT^2 + dT^3 + eT^4;

    let density_value_kg_per_m3 = 
        a 
        + b * fluid_temp_kelvin;


    return Ok(MassDensity::new::<
              kilogram_per_cubic_meter>(density_value_kg_per_m3));
}

/// function to obtain nitrate salt viscosity
/// given a temperature
///
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68./// Jana, S. S., 
/// Maheshwari, N. K., & Vijayan, P. K. (2016). 
///
/// mu Pa-s (T = 440 - 500 K) 
/// = 0.93845
/// -0.54754 T(K)
/// + 1.08225e-5 T(K)^2
/// - 7.2058e-9 T(K)^3
///
/// mu Pa-s (T = 500 - 800 K) 
/// = 0.23816
/// - 1.2768e-3 T(K)
/// + 2.6275e-6 T(K)^2
/// - 2.4331e-9 T(K)^3
/// + 8.507e-13 T(K)^4
///
/// 
pub fn get_hitec_equimolar_viscosity(
    fluid_temp: ThermodynamicTemperature) -> Result<DynamicViscosity,
ThermalHydraulicsLibError>{

    range_check_hitec_salt(fluid_temp)?;
    let fluid_temp_kelvin = fluid_temp.get::<kelvin>();
    // mu Pa-s (T = 500 - 800 K) 
    // = 0.23816
    // - 1.2768e-3 T(K)
    // + 2.6275e-6 T(K)^2
    // - 2.4331e-9 T(K)^3
    // + 8.507e-13 T(K)^4
    let mut a = 0.23816;
    let mut b = - 1.2768e-3;
    let mut c = 2.6275e-6;
    let mut d = -2.4331e-9;
    let mut e = 8.507e-13;

    if fluid_temp_kelvin < 500.0 {
        // mu Pa-s (T = 440 - 500 K) 
        // = 0.93845
        // -0.54754 T(K)
        // + 1.08225e-5 T(K)^2
        // - 7.2058e-9 T(K)^3
        a = 0.93845;
        b = - 5.4754e-3;
        c = 1.08225e-5;
        d = -7.2058e-9;
        e = 0.0;

    }

    // generic correlation is:
    // a + bT + cT^2 + dT^3 + eT^4;
    let viscosity_value_pascal_second = 
        a 
        + b * fluid_temp_kelvin
        + c * fluid_temp_kelvin.powf(2.0)
        + d * fluid_temp_kelvin.powf(3.0)
        + e * fluid_temp_kelvin.powf(4.0);

    Ok(DynamicViscosity::new::<pascal_second>(viscosity_value_pascal_second))
                                
}

/// function to obtain nitrate salt specific heat capacity
/// given a temperature
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68./// Jana, S. S., 
/// Maheshwari, N. K., & Vijayan, P. K. (2016). 
///
/// cp (J/kg/K) = 1443.0 + 0.172 T
/// T in degc
pub fn get_hitec_constant_pressure_specific_heat_capacity(
    fluid_temp: ThermodynamicTemperature) -> Result<SpecificHeatCapacity,
ThermalHydraulicsLibError>{

    range_check_hitec_salt(fluid_temp)?;
    // note, specific entropy and heat capcity are the same unit...
    //
    let temperature_degrees_c_value = fluid_temp.get::<degree_celsius>();
    let cp_value_joule_per_kg = 1443.0 + 0.172*temperature_degrees_c_value;

    Ok(SpecificHeatCapacity::new::<joule_per_kilogram_kelvin>(
        cp_value_joule_per_kg))
}

/// function to obtain nitrate salt thermal conductivity
/// given a temperature
///
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68./// Jana, S. S., 
/// Maheshwari, N. K., & Vijayan, P. K. (2016). 
///
/// k (thermal conductivity in W/mK) = 0.443 + 1.9e-4 T
/// T in degc
pub fn get_hitec_thermal_conductivity(
    fluid_temp: ThermodynamicTemperature) -> Result<ThermalConductivity,ThermalHydraulicsLibError> {


    range_check_hitec_salt(fluid_temp)?;
    let thermal_conductivity_value = 0.443 - 1.9e-4* fluid_temp
        .get::<degree_celsius>();

    return Ok(ThermalConductivity::new::<watt_per_meter_kelvin>(
        thermal_conductivity_value));
}


/// function checks if a fluid temperature falls in a range 
///
/// If it falls outside this range, it will panic
/// or throw an error, and the program will not run
///
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68./// Jana, S. S., 
/// Maheshwari, N. K., & Vijayan, P. K. (2016). 
///
/// From HITEC, the applicable range is 440K - 800 K, 
///
/// In Du's paper, the viscosity correlation is applicable from 440 to 800K
/// while the rest of the properties are from 420-800K
/// 
///
pub fn range_check_hitec_salt(fluid_temp: ThermodynamicTemperature) 
    -> Result<bool,ThermalHydraulicsLibError>{

        // first i convert the fluidTemp object into a degree 
        // celsius

        range_check(&Material::Liquid(LiquidMaterial::DowthermA), 
            fluid_temp, 
            ThermodynamicTemperature::new::<kelvin>(800.0), 
            ThermodynamicTemperature::new::<kelvin>(440.0))?;

        return Ok(true);

    }
