use ndarray::*;
use ndarray_linalg::error::LinalgError;
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
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

/// this is essentially a 1D pipe array containing two CVs 
/// and two other radially connected arrays
///
/// Usually, these will be nested inside a heat transfer component 
/// and then be used
///
/// Within this array, the implicit Euler Scheme is used
#[derive(Debug,Clone,PartialEq,Default)]
pub struct FluidArray {

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

    /// number of inner nodes in the array besides the front and back 
    /// cv 
    ///
    /// so total number of nodes is at least 2 all the time 
    /// as the front and back CV count as nodes 
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

}

impl FluidArray {


    /// obtains a clone of the temperature vector within the CV 
    /// thus obtaining the temperature profile
    pub fn get_temperature_vector(&self) -> Result<
    Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        let mut temperature_vec: Vec<ThermodynamicTemperature> = vec![];

        for temperature in self.temperature_array_current_timestep.iter() {
            temperature_vec.push(temperature.clone());
        }

        return Ok(temperature_vec);
    }

    /// sets the temperature vector to a 
    pub fn set_temperature_vector(&mut self,
    temperature_vec: Vec<ThermodynamicTemperature>) -> Result<(), ThermalHydraulicsLibError>{

        let number_of_temperature_nodes = self.inner_nodes + 2;

        // check if temperature_vec has the correct number_of_temperature_nodes

        if temperature_vec.len() !=  number_of_temperature_nodes {
            let shape_error = ShapeError::from_kind(
                ErrorKind::IncompatibleShape
            );

            let linalg_error = LinalgError::Shape(shape_error);

            return Err(ThermalHydraulicsLibError::LinalgError
                (linalg_error));

        }

        let mut temperature_array: Array1<ThermodynamicTemperature>
        = Array::default(number_of_temperature_nodes);

        for (index,temperature) in temperature_array.iter_mut().enumerate() {
            *temperature = temperature_vec[index];
        }
        
        Ok(())
    }
}
