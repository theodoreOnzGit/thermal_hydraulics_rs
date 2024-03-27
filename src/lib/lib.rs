// Note: //! indicates crate level documentation
//
//
//! A Library which contains useful traits and methods for thermal 
//! hydraulics calculations.
//!
//!
//! This crate has heavy reliance on units of measure (uom) released under 
//! Apache 2.0 license. So you'll need to get used to unit safe calculations
//! with uom as well.
//!
//!
//! This library was initially developed for 
//! use in my PhD thesis under supervision 
//! of Professor Per F. Peterson. It a thermal hydraulics
//! library in Rust that is released under the GNU General Public License
//! v 3.0. This is partly due to the fact that some of the libraries 
//! inherit from GeN-Foam and OpenFOAM, both licensed under GNU General
//! Public License v3.0.
//!
//! As such, the entire library is released under GNU GPL v3.0. It is a strong 
//! copyleft license which means you cannot use it in proprietary software.
//!
//!
//! License
//!    This is a thermal hydraulics library written 
//!    in rust meant to help with the
//!    fluid mechanics and heat transfer aspects of the calculations
//!    for the Compact Integral Effects Tests (CIET) and hopefully 
//!    Gen IV Reactors such as the Fluoride Salt cooled High Temperature 
//!    Reactor (FHR)
//!     
//!    Copyright (C) 2022-2023  Theodore Kay Chen Ong, Singapore Nuclear
//!    Research and Safety Initiative, Per F. Peterson, University of 
//!    California, Berkeley Thermal Hydraulics Laboratory
//!
//!    thermal_hydraulics_rs is free software; you can 
//!    redistribute it and/or modify it
//!    under the terms of the GNU General Public License as published by the
//!    Free Software Foundation; either version 2 of the License, or (at your
//!    option) any later version.
//!
//!    thermal_hydraulics_rs is distributed in the hope 
//!    that it will be useful, but WITHOUT
//!    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//!    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//!    for more details.
//!
//!    This thermal hydraulics library 
//!    contains some code copied from GeN-Foam, and OpenFOAM derivative.
//!    This offering is not approved or endorsed by the OpenFOAM Foundation nor
//!    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
//!    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.
//!    Nor is it endorsed by the authors and owners of GeN-Foam.
//!
//!    You should have received a copy of the GNU General Public License
//!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//!
//! © All rights reserved. Theodore Kay Chen Ong,
//! Singapore Nuclear Research and Safety Initiative,
//! Per F. Peterson,
//! University of California, Berkeley Thermal Hydraulics Laboratory
//!
//! Main author of the code: Theodore Kay Chen Ong, supervised by
//! Professor Per F. Peterson
//!
//! Btw, I no affiliation with the Rust Foundation. 
//!
#![warn(missing_docs)]
extern crate uom;

/// Fluid Mechanics Module (testing)
pub mod fluid_mechanics_lib;

/// Heat Transfer Module (testing)
pub mod heat_transfer_lib;

/// for mostly incompressible fluids using the Boussinesq Approximation
/// that is, density doesn't change much except for natural convection
pub mod boussinesq_solver;

/// use peroxide macros 
#[macro_use]
extern crate peroxide;

/// provides error types for thermal_hydraulics_rs
pub mod thermal_hydraulics_error;

/// prelude, for easy importing 
pub mod prelude;

// to do:
// 1. transfer heat transfer sandbox to use boussinesq solver 
// 2. test the infinite medium test case using the pre built components
// 3. build and validate natural circulation loop
// 4. update the steel ss304L libraries to use correlation rather than spline building
// this is because graves correlation steel temperature range is too small
