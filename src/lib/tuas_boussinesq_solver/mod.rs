//#![warn(missing_docs)]
///// Module specifically for thermophysical properties
///// For liquids and solids with almost invariable density
/////
//pub mod boussinesq_thermophysical_properties;
//
///// Module for correlations of fluid mechanics 
///// suitable for tuas_boussinesq_solver (single phase, negligble density changes
///// except for buoyancy)
//pub mod fluid_mechanics_correlations;
//
///// Module for heat transfer correlations 
///// suitable for tuas_boussinesq_solver (single phase, negligble density changes
///// except for buoyancy)
//pub mod heat_transfer_correlations;
//
///// specific dimensions for control volume construction
//pub mod control_volume_dimensions;
//
///// Module for boundary conditions 
//pub mod boundary_conditions;
//
///// Module for single control volumes (mainly for fluid control volumes,
///// but solid control volumes are set by setting flowrate to zero)
/////
///// Single control volumes by default have functions which abstract away 
///// the details of calculating heat transfer between different 
///// single control volumes as well as between single control volumes and 
///// different boundary conditions
/////
///// This will abstract away some functionality of the following 
///// modules, and is therefore dependent on these modules:
/////
///// 1. boussinesq_thermophysical_properties
///// 2. fluid_mechanics_correlations
///// 3. heat_transfer_correlations
///// 4. control_volume_dimensions
///// 5. boundary_conditions
/////
///// By itself, it will NOT contain functions on how to interact with array 
///// control volumes. This is to prevent overbloated hard to read code
//pub mod single_control_vol;
//
//
///// Module for array control volumes (mainly for fluid control volumes,
///// but solid control volumes are set by setting flowrate to zero)
///// suitable for tuas_boussinesq_solver (single phase, negligble density changes
///// except for buoyancy)
/////
///// also contains code to help calculate pressure drop and mass flow rate 
///// amongst multiple fluid components (eg. pipes) which are usually 
///// represented by array control volumes
///// This will abstract away some functionality of the following 
///// modules, and is therefore dependent on these modules:
/////
///// 1. boussinesq_thermophysical_properties
///// 2. fluid_mechanics_correlations
///// 3. heat_transfer_correlations
///// 4. control_volume_dimensions
///// 5. boundary_conditions
///// 6. single_control_vol
/////
///// By itself, it will NOT contain functions on how to interact with array 
///// control volumes. This is to prevent overbloated hard to read code
//pub mod array_control_vol_and_fluid_component_collections;
//
///// Module for pre-built-components 
///// suitable for tuas_boussinesq_solver (single phase, negligble density changes
///// except for buoyancy)
/////
///// It's dependent on all the other modules within the tuas_boussinesq_solver
/////
///// You don't want to write everything from scratch right? 
//pub mod pre_built_components;
