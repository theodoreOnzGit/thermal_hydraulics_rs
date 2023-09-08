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


/// This module contains direct translations of GeN-Foam code 
/// into rust for the lumped nuclear structure
///
/// this is a useful reference to start from for finite difference or 
/// finite volume implicitly solved code 
#[path = "./gen-foam-lumped-nuclear-structure.rs"]
pub mod gen_foam_lumped_nuclear_structure;
/// This module contains adapts the GeN-Foam code for the 
/// lumped nuclear structure
/// for the heat transfer module
pub mod lumped_nuclear_structure_inspired_functions;

/// this module contains constructors for an array cv which 
/// functions as a one dimensional cartesian conduction medium
pub mod one_dimension_cartesian_conducting_medium;
pub use one_dimension_cartesian_conducting_medium::*;



/// contains functions to advance timestep and such
pub mod calculation;


/// contains functions to get max timestep and other things 
/// also deals with preparatory steps before calculation 
/// for example, linking up control volumes 
///
/// linking of control volumes
pub mod preprocessing;
pub use preprocessing::*;

/// contains functions or methods to get temperature 
/// and other things 
pub mod postprocessing;
pub use postprocessing::*;



/// sandbox, for miscellaneous testing of code
mod sandbox;

/// contains matrix calculations specific to fluid nodes
/// arranged in a 1D array
pub mod fluid_nodes;

/// contains matrix calculations specific to solid nodes 
/// these are meant to represent the "shell" of the pipe 
/// or any kind of solid material in the pipe
pub mod solid_nodes;

/// one dimensional pipe 
/// it's a control volume for one dimensional pipes 
/// will describe ONLY the fluid parts without the solid shells 
pub mod one_dimension_fluid_array;

/// one dimensional solid rod like object 
/// we call it a column to generalise the idea that is should 
/// be somewhat cylindrical, but may include other shapes
/// However, the cross sectional area is uniform
pub mod one_dimension_solid_array;
