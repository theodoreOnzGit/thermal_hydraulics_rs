/// for solid-solid and solid-fluid interaction, heat transfer can 
/// be expressed in terms of thermal conductance or thermal 
/// resistance 
///
/// this module contains functions for these 
pub mod conductance;
pub use conductance::*;

/// for fluid-fluid interaction, where fluid flows and carries 
/// heat from one control volume to another, we call this advection 
/// functions to calculate advection and placed in this module
pub mod advection;
pub use advection::*;

