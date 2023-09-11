
/// calculates an appropriate time step for two single control volume 
/// nodes, appending the max time step allowable for each one
/// first based on thermal conductivity
///
/// also should have courant number (TBC)
pub mod two_single_cv;
pub use two_single_cv::*;



/// Calculates a time step based on fourier number for a single 
/// control volume, usually tied to some boundary condition
///
pub mod cv_and_bc;
pub use cv_and_bc::*;
