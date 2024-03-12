#![warn(missing_docs)]
/// Module specifically for thermophysical properties
/// For liquids and solids with almost invariable density
///
pub mod boussinesq_thermophysical_properties;

/// Module for correlations of fluid mechanics 
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod fluid_mechanics_correlations;


/// Module for heat transfer correlations 
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod heat_transfer_correlations;

/// Module for solid control volumes (both singular control volumes 
/// and arrays of control volumes)
pub mod solid_control_vol;

/// Module for fluid control volumes (both singular control volumes 
/// and arrays of control volumes)
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod fluid_control_vol;

/// Module for pre-built-components 
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod pre_built_components;
