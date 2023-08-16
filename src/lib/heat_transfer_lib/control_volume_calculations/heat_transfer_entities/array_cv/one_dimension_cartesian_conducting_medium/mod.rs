use ndarray::*;
use uom::si::f64::*;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::HeatTransferEntity;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::
thermophysical_properties::Material;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::CVType;

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
}

/// here, i mostly do constructors
impl CartesianConduction1DArray {

    /// constructs a new instance of the CartesianConduction1DArray
    pub fn new(material: Material,
    initial_uniform_temperature: ThermodynamicTemperature,
    uniform_pressure: Pressure,
    inner_nodes: usize,
    total_length: Length) -> HeatTransferEntity {
        // we start building the 1Darray object by a default first
        let mut array_to_return = Self::default();

        // set the scalars first, they are the easiest
        array_to_return.material_control_volume = material;
        array_to_return.inner_nodes = inner_nodes;
        array_to_return.pressure_control_volume = uniform_pressure;
        array_to_return.total_length = total_length;

        // by now, we should be able to set the volume fraction array 
        let vol_frac_array = 
            array_to_return.construct_volume_fraction_array().unwrap();

        array_to_return.volume_fraction_array = vol_frac_array.clone();

        // set the temperature arrays
        let number_of_temperature_nodes: usize = inner_nodes + 2;

        let mut initial_temperature_array: 
        Array1<ThermodynamicTemperature> = 
            Array::default(number_of_temperature_nodes);

        initial_temperature_array.fill(initial_uniform_temperature);

        array_to_return.temperature_array_current_timestep = 
            initial_temperature_array.clone();

        array_to_return.temperature_array_next_timestep = 
            initial_temperature_array;

        // lets set the remaining control volumes
        // 
        // for a 1D volume, the volume fraction is the same 
        // as the length fraction
        let boundary_length: Length = vol_frac_array[0] * total_length;

        let boundary_cv_entity: HeatTransferEntity = 
        SingleCVNode::new_one_dimension_volume(
            boundary_length,
            material,
            initial_uniform_temperature,
            uniform_pressure).unwrap();

        let boundary_cv: SingleCVNode = match boundary_cv_entity {
            HeatTransferEntity::ControlVolume(cv) => {

                match cv {
                    CVType::SingleCV(cv) => {
                        cv
                    },
                    _ => panic!(),
                }

            },
            _ => panic!(),
        };

        array_to_return.inner_single_cv = boundary_cv.clone();
        array_to_return.outer_single_cv = boundary_cv.clone();

        // now package as an array cv 

        let heat_transfer_entity 
        = HeatTransferEntity::ControlVolume(
            CVType::ArrayCV(
                ArrayCVType::Cartesian1D(array_to_return)
            )
        );

        return heat_transfer_entity;
    }



}

/// Functions or methods to retrieve temperature and other such 
/// data from the array_cv
pub mod postprocessing;
pub use postprocessing::*;

/// Functions or methods to get timestep and other such quantiies 
/// for calculations 
///
/// helps to set up quantities used in calculation step
pub mod preprocessing;
pub use preprocessing::*; 

/// Contains functions which advance the timestep
/// it's the bulk of calculation
pub mod calculation;
pub use calculation::*;
