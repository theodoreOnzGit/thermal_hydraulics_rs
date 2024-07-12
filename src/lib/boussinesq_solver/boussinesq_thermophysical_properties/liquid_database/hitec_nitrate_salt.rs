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
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::dynamic_viscosity::millipascal_second;
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
/// rho (kg/m3) = 2090 - 0.636 T(C)
pub fn get_hitec_equimolar_density(
    fluid_temp: ThermodynamicTemperature) -> Result<MassDensity,ThermalHydraulicsLibError> {


    // first we check if fluid temp is between 220-630C (range of validity)
    // panic otherwise
    range_check_hitec_salt(fluid_temp)?;

    //then convert the fluidTemp object into a f64
    // and plug it into the correlation
    let density_value_kg_per_m3 = 2090.0 - 0.636*fluid_temp
       .get::<degree_celsius>();

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
/// mu mPa-s = 22.714 - 0.120 T + 2.281e-4 T^2 -1.474e-7 T^3
/// T in degc
pub fn get_hitec_equimolar_viscosity(
    fluid_temp: ThermodynamicTemperature) -> Result<DynamicViscosity,
ThermalHydraulicsLibError>{

    range_check_hitec_salt(fluid_temp)?;
    let temperature_degrees_c_value = fluid_temp.get::<degree_celsius>();
    let viscosity_value_millipascal_second = 
        22.714 
        - 0.120*temperature_degrees_c_value.powf(1.0)
        + 2.281e-4*temperature_degrees_c_value.powf(2.0)
        - 1.474e-7*temperature_degrees_c_value.powf(3.0);

    Ok(DynamicViscosity::new::<millipascal_second>(viscosity_value_millipascal_second))
                                
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
/// From HITEC, the applicable range is 440C - 800 C, 
///
/// In Du's paper, the viscosity correlation is applicable from 440 to 800C
/// while the rest of the properties are from 420-800C
/// 
///
pub fn range_check_hitec_salt(fluid_temp: ThermodynamicTemperature) 
    -> Result<bool,ThermalHydraulicsLibError>{

        // first i convert the fluidTemp object into a degree 
        // celsius

        range_check(&Material::Liquid(LiquidMaterial::DowthermA), 
            fluid_temp, 
            ThermodynamicTemperature::new::<degree_celsius>(800.0), 
            ThermodynamicTemperature::new::<degree_celsius>(440.0))?;

        return Ok(true);

    }
