/// suppose the control volume interacts with a BC which is 
/// a constant heat addition, which is a constant power rating
/// 
///
/// cooling is also possible, just supply a negative Power 
/// quantity
/// 
/// for the most part, this is a wrapper function
pub mod constant_heat_addition;
pub use constant_heat_addition::*;


/// suppose control volume interacts with a constant heat flux BC
/// we will need a sort of surface area in order to determine the 
/// power added
pub mod constant_heat_flux;
pub use constant_heat_flux::*;

/// suppose control volume interacts with constant temperature BC
/// we need some thermal conductance value to obtain a power
/// value
pub mod constant_temperature;
pub use constant_temperature::*;

// this function specifically calculates interaction 
// between two single CV nodes
// 
pub mod two_single_cv;
pub use two_single_cv::*;
