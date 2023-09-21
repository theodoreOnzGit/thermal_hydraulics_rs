//! this module contains functions to calculate the enthalpy
//! of a control volume at the the next timestep
//!
//!
//! For Control Volume calculations in general
//! we have the form:
//!
//! dH_cv/dt = H_in - H_out + Q_s + W_s
//!
//! H_cv is the control volume enthalpy
//!
//! H_in is the sum of enthalpy flows in
//!
//! H_out is the sum of enthalpy flows out
//!
//! Q_s is the heat supplied to the volume per second
//!
//! W_s is the work done on the system per second
//!
//!
//! After discretisation, we can use:
//!
//! (H_cv (t+1) - H_cv (t)) / dt = H_in - H_out + Q_s + W_s
//!
//! H_cv (t+1)  = dt * (H_in - H_out + Q_s + W_s) + H_cv (t)
//!
//! It remains to be seen whether the enthalpy flows in and
//! out are calculated at the current time step  (explicit)
//! or next time step (implict)
//!
//! Of course, implicit calculations are more stable but
//! slower in general than explicit calculations
//!
//! we will be using the uom module to ensure that calculations are
//! done with correct units

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






/// This module contains commonly used functions for Explicit and
/// Implicit timestep calculations
pub mod common_functions;

/// heat_transfer_entities contain the structs and enums which describe 
/// control volumes and boundary conditions
pub mod heat_transfer_entities;

/// heat_transfer_interactions contains the structs and enums 
/// which describe behaviour for how control volumes and BCs 
/// (heat_transfer_entities) interact with each other
pub mod heat_transfer_interactions;

///// This module contains traits useful for constructing control volumes
/////
///// The difference between this and the fluid entity structs and traits
///// is that the control volume traits are much more generic,
///// there can be any number of fluids flowing in and out of the control 
///// volume whereas FluidEntity_StructsAndTraits deals with single
///// flow input and output coming out of the pipe
/////
///// Furthermore, structs and traits found in FluidEntity_StructsAndTraits
///// tries to incorporate the steps required to solve for temperature at
///// every time step. Control volume traits do not help you in this manner
/////
//mod control_volume_traits;

///// This module contains Structs and Traits for
///// Generic Fluid entities
/////
///// It is still work in progress, so not quite done yet
/////
/////
/////
//#[allow(non_snake_case)]
//mod FluidEntity_StructsAndTraits;

///// This module just helps document the development process and sandboxing
///// i have done
///// it serves more as a scratchpad
///// stuff in here should not be inherited
//#[allow(non_snake_case)]
//mod Sandbox;



///// This module contains functions which help to calculate
///// the enthalpy explicitly, ie using enthalpy in and out for current
///// timestep
/////
///// The known information from which to start is 
///// (1) the mass flow rate
///// (2) the temperatures of each part of the fluid
///// (3) mass of the control volume
//
//#[allow(non_snake_case)]
//mod ExplicitCalculations;

///// This module contains implementations for therminol VP1
///// pipe heat transfer
//#[allow(non_snake_case)]
//mod TherminolDowthermPipes;





