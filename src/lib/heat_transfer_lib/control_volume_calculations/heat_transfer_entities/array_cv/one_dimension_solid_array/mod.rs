use ndarray::*;
use uom::si::f64::*;

use crate::fluid_mechanics_lib::fluid_component_calculation::enums::DimensionlessDarcyLossCorrelations;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::nusselt_correlations::enums::NusseltCorrelation;
use crate::heat_transfer_lib::
thermophysical_properties::Material;

/// this is essentially a 1D pipe array containing two CVs 
/// and two other laterally connected arrays
/// (it's essentially a generic solid array representing heat 
/// structures with mainly axial conduction and radial conduction)
///
/// it can be used to represent rods, or cylindrical shells
/// in the latter case, the Column is hollow so to speak
///
/// Usually, these will be nested inside a heat transfer component 
/// and then be used
///
/// Within this array, the implicit Euler Scheme is used
///
/// You must supply the number of nodes for the fluid array
/// Note that the front and back cv count as one node
#[derive(Debug,Clone,PartialEq)]
pub struct SolidColumn {

    /// represents the control volume at the back 
    /// imagine a car or train cruising along in a positive x direction
    ///
    /// //----------------------------------------------> x 
    ///
    /// //            (back --- train/car --- front)
    /// //            lower x                 higher x
    ///
    pub back_single_cv: SingleCVNode,

    /// represents the control volume at the front
    ///
    /// to think of which is front and back, we think of coordinates 
    /// imagine a car or train cruising along in a positive x direction
    ///
    /// //----------------------------------------------> x 
    ///
    /// //            (back --- train/car --- front)
    /// //            lower x                 higher x
    ///
    pub front_single_cv: SingleCVNode,

    /// number of inner nodes ,
    /// besides the back and front node 
    ///
    /// total number of nodes is inner_nodes + 2
    inner_nodes: usize,


    // total length for the array
    total_length: Length,

    // cross sectional area for the 1D array, assumed to be uniform 
    xs_area: Area,

    /// temperature array current timestep 
    /// only accessible via get and set methods
    temperature_array_current_timestep: Array1<ThermodynamicTemperature>,

    // temperature_array_next timestep 
    temperature_array_next_timestep: Array1<ThermodynamicTemperature>,

    /// control volume material 
    pub material_control_volume: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,

    // volume fraction array 
    volume_fraction_array: Array1<f64>,


    /// now solid arrays (columns) can be connected to solid arrays 
    /// or other fluid arrays adjacent to it radially
    ///
    /// There will be no advection but there can 
    /// be thermal conductance shared between the nodes 
    ///
    /// hence, I only want to have a copy of the temperature 
    /// arrays radially adjacent to it
    ///
    /// plus their thermal resistances
    /// N is the array size, which is known at compile time

    pub lateral_adjacent_array_temperature_vector: 
    Vec<Array1<ThermodynamicTemperature>>,

    /// now solid arrays (columns) can be connected to solid arrays 
    /// or other fluid arrays adjacent to it radially
    ///
    /// There will be no advection 
    /// but there can be thermal conductance shared between the nodes 
    ///
    /// hence, I only want to have a copy of the temperature 
    /// arrays radially adjacent to it
    ///
    /// plus their thermal resistances 
    /// N is the array size, which is known at compile time
    pub lateral_adjacent_array_conductance_vector:
    Vec<Array1<ThermalConductance>>,

    /// solid arrays can also be connected to heat sources 
    /// or have specified volumetric heat sources 
    pub q_vector: Vec<Power>,

    /// solid arrays should have their power distributed according 
    /// to their nodes 
    pub q_fraction_vector: Vec<Array1<f64>>,

}
