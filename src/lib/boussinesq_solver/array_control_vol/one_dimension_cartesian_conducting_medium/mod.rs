use ndarray::*;
use uom::si::f64::*;

use crate::boussinesq_solver::single_control_vol::SingleCVNode;
use crate::boussinesq_solver::heat_transfer_correlations::nusselt_number_correlations::enums::NusseltCorrelation;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::Material;


/// for 1D Cartesian Conduction array,
/// it is essentially an array control volume of one homogeneous 
/// material
///
/// it is in Cartesian coordinates, basically, x direction only conduction 
///
/// the structure is segregated into several smaller nodes using finite 
/// difference methods 
///
/// I'll use lumped_nuclear_structure_inspired_functions to calculate 
/// new temperatures for this structure 
///
/// the scheme used is the implicit Euler scheme to calculate new 
/// temperatures. However, material properties are calculated using 
/// current timestep temperatures rather than next timestep temperatures 
/// therefore, it is more of a hybrid between the implicit and explicit 
/// schemes. 
///
/// the important methods are to advance timestep, and to update 
/// material properties at every timestep
#[derive(Debug,Clone,PartialEq,Default)]
pub struct CartesianConduction1DArray {


    /// represents the inner (lower r) control volume end 
    /// or back (lower x) control volume
    /// to think of which is front and back, we think of coordinates 
    /// imagine a car or train cruising along in a positive x direction
    ///
    /// //----------------------------------------------> x 
    ///
    /// //            (back --- train/car --- front)
    /// //            lower x                 higher x
    ///
    pub inner_single_cv: SingleCVNode,
    /// represents the outer (higher r) control volume end 
    /// or front (higher x) control volume
    ///
    /// to think of which is front and back, we think of coordinates 
    /// imagine a car or train cruising along in a positive x direction
    ///
    /// //----------------------------------------------> x 
    ///
    /// //            (back --- train/car --- front)
    /// //            lower x                 higher x
    ///
    pub outer_single_cv: SingleCVNode,

    // number of ADDITIONAL nodes within the array 
    // in addition to the inner and outer single cvs
    inner_nodes: usize,
    
    // total length for the 1D array
    total_length: Length,

    // volume fraction array 
    volume_fraction_array: Array1<f64>,

    /// temperature array current timestep 
    pub temperature_array_current_timestep: Array1<ThermodynamicTemperature>,

    // temperature_array_next timestep 
    temperature_array_next_timestep: Array1<ThermodynamicTemperature>,

    /// control volume material 
    pub material_control_volume: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,

    /// nusselt correlation type 
    nusselt_correlation_type: NusseltCorrelation
}
