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



//extern crate uom;
//use uom::si::f64::*;




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





