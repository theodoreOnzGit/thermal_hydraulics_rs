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
use uom::si::dynamic_viscosity::pascal_second;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::available_energy::joule_per_kilogram;

// this is for the root finding algorithms
extern crate peroxide;
use peroxide::prelude::*;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::{range_check, LiquidMaterial, Material};
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

pub fn get_flinak_density(
    fluid_temp: ThermodynamicTemperature) -> Result<MassDensity,ThermalHydraulicsLibError> {
    todo!()
}

pub fn get_flinak_dynamic_viscosity(
    fluid_temp: ThermodynamicTemperature) -> Result<DynamicViscosity,
ThermalHydraulicsLibError>{
    todo!()
}

pub fn get_flinak_constant_pressure_specific_heat_capacity(
    fluid_temp: ThermodynamicTemperature) -> Result<SpecificHeatCapacity,
ThermalHydraulicsLibError>{
    todo!()
}

pub fn get_flinak_thermal_conductivity(
    fluid_temp: ThermodynamicTemperature) -> Result<ThermalConductivity,ThermalHydraulicsLibError> {
    todo!()
}
pub fn get_flinak_specific_enthalpy(
    fluid_temp: ThermodynamicTemperature) -> 
Result<AvailableEnergy,ThermalHydraulicsLibError>{
    todo!()
}
pub fn get_temperature_from_enthalpy(
    fluid_enthalpy: AvailableEnergy) -> Result<ThermodynamicTemperature,ThermalHydraulicsLibError> {

    // if enthalpy value below zero,
    // based on me setting zero enthalpy at the lower end of the 
    // temperature validity range for enthalpy,
    // then enthalpy is technically out of range
    if fluid_enthalpy.value < 0_f64 {
        panic!("FLiNaK : get_temperature_from_enthalpy \n
               enthalpy < 0.0 , out of correlation range");
    }

    // first let's convert enthalpy to a double (f64)
    let enthalpy_value_joule_per_kg = 
        fluid_enthalpy.get::<joule_per_kilogram>();

    // second let's define a function 
    // or actually a closure or anonymous function that
    // is aware of the variables declared
    // LHS is actual enthalpy value

    let enthalpy_root = |temp_degrees_kelvin_value : AD| -> AD {
        let lhs_value = enthalpy_value_joule_per_kg;
        // convert AD type into double
        let temp_degrees_kelvin_value_double = temp_degrees_kelvin_value.x();

        let fluid_temperature = 
            ThermodynamicTemperature::new::<kelvin>(
                temp_degrees_kelvin_value_double);
        let rhs = get_flinak_specific_enthalpy(fluid_temperature).unwrap();
        let rhs_value = rhs.get::<joule_per_kilogram>();

        return AD0(lhs_value-rhs_value);
    };
    
    // now solve using bisection
    // the range is from 732.2 K - 1573 K
    
    let fluid_temperature_degrees_kelvin_result 
        = bisection(enthalpy_root,
                    (732.2,1573.0),
                    100,
                    1e-8);

    let fluid_temperature_degrees_kelvin = fluid_temperature_degrees_kelvin_result.unwrap();

    return Ok(ThermodynamicTemperature::
        new::<kelvin>(fluid_temperature_degrees_kelvin));

}
/// function checks if a fluid temperature falls in a range 
///
/// If it falls outside this range, it will panic
/// or throw an error, and the program will not run
///
/// Sohal, M. S., Ebner, M. A., Sabharwall, P., & Sharpe, P. (2010). 
/// Engineering database of liquid salt thermophysical and thermochemical 
/// properties (No. INL/EXT-10-18297). Idaho National Lab.(INL), 
/// Idaho Falls, ID (United States).
///
/// For FLiBe, the applicable range is 732.2K (melting point) - 1573 K.
/// I try to make the range as wide as possible because Gnielinski's correlation 
/// requires corrections using wall temperature. These may be outside 
/// the usual bulk temperatures of FLiBe.
///
/// 
pub fn range_check_flinak_salt(fluid_temp: ThermodynamicTemperature) 
    -> Result<bool,ThermalHydraulicsLibError>{

        // first i convert the fluidTemp object into a degree 
        // celsius

        range_check(&Material::Liquid(LiquidMaterial::FLiBe), 
            fluid_temp, 
            ThermodynamicTemperature::new::<kelvin>(1573.0), 
            ThermodynamicTemperature::new::<kelvin>(732.2))?;

        return Ok(true);

    }

