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

/// Module for boundary conditions 
pub mod boundary_conditions;

/// Module for single control volumes (mainly for fluid control volumes,
/// but solid control volumes are set by setting flowrate to zero)
///
/// Single control volumes by default have functions which abstract away 
/// the details of calculating heat transfer between different 
/// single control volumes as well as between single control volumes and 
/// different boundary conditions
///
/// This is, it will abstract away some functionality of the following 
/// modules, and is therefore dependent on these modules:
///
/// 1. boussinesq_thermophysical_properties
/// 2. fluid_mechanics_correlations
/// 3. heat_transfer_correlations
/// 4. control_volume_dimensions
/// 5. boundary_conditions
///
/// By itself, it will NOT contain functions on how to interact with array 
/// control volumes. This is to prevent overbloated hard to read code
pub mod single_control_vol;


/// Module for array control volumes (mainly for fluid control volumes,
/// but solid control volumes are set by setting flowrate to zero)
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod array_control_vol;

/// Module for pre-built-components 
/// suitable for boussinesq_solver (single phase, negligble density changes
/// except for buoyancy)
pub mod pre_built_components;
