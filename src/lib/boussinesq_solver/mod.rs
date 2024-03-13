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

/// specific dimensions for control volume construction
pub mod control_volume_dimensions;

/// Module for single control volumes (mainly for fluid control volumes,
/// but solid control volumes are set by setting flowrate to zero)
pub mod single_control_vol;

/// Module for boundary conditions 
pub mod boundary_conditions;

/// Module for array control volumes (mainly for fluid control volumes,
/// but solid control volumes are set by setting flowrate to zero)
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod array_control_vol;

/// Module for pre-built-components 
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod pre_built_components;
