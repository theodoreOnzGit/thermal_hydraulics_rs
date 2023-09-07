use ndarray::*;
use ndarray_linalg::error::LinalgError;
use uom::si::f64::*;

use crate::fluid_mechanics_lib::fluid_component_calculation::enums::DimensionlessDarcyLossCorrelations;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::HeatTransferEntity;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::nusselt_correlations::enums::NusseltCorrelation;
use crate::heat_transfer_lib::
thermophysical_properties::Material;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::CVType;
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::specific_enthalpy;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

/// this is essentially a 1D pipe array containing two CVs 
/// and two other radially connected arrays
/// (it's essentially a generic array representing pipes or 
/// other fluid components with one inlet and one outlet)
///
/// Usually, these will be nested inside a heat transfer component 
/// and then be used
///
/// Within this array, the implicit Euler Scheme is used
///
/// You must supply the number of nodes for the fluid array
/// it will then be able to construct fixed size temperature 
/// arrays. Note that the front and back cv count as one node
///
#[derive(Debug,Clone,PartialEq)]
pub struct FluidArray<const NUMBER_OF_NODES: usize> {

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


    // total length for the array
    total_length: Length,

    // cross sectional area for the 1D array, assumed to be uniform 
    xs_area: Area,

    /// temperature array current timestep 
    /// only accessible via get and set methods
    temperature_array_current_timestep: [ThermodynamicTemperature; NUMBER_OF_NODES],

    // temperature_array_next timestep 
    temperature_array_next_timestep: [ThermodynamicTemperature; NUMBER_OF_NODES],

    /// control volume material 
    pub material_control_volume: Material,

    /// for fluid nodes, it is surrounded by some solid material,
    /// store it here 
    pub adjacent_solid_material: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,

    // volume fraction array 
    volume_fraction_array: [f64; NUMBER_OF_NODES],

    /// mass flowrate through the fluid array, 
    mass_flowrate: MassRate,

    /// pressure loss term 
    pressure_loss: Pressure,

    /// wetted perimeter (for hydraulic diameter) 
    wetted_perimiter: Length,

    /// incline angle 
    incline_angle: Angle,

    /// internal pressure source 
    internal_pressure_source: Pressure,

    /// pipe loss properties 
    pipe_loss_properties: DimensionlessDarcyLossCorrelations,

    /// nusselt correlation 
    nusselt_correlation: NusseltCorrelation,

    /// now fluid arrays can be connected to solid arrays 
    /// or other fluid arrays adjacent to it radially
    ///
    /// There will be no advection in the radial direction,
    /// but there can be thermal conductance shared between the nodes 
    ///
    /// hence, I only want to have a copy of the temperature 
    /// arrays radially adjacent to it
    ///
    /// plus their thermal resistances (one average will be given 
    /// per array for expedience)
    /// N is the array size, which is known at compile time

    radial_adjacent_array_temperature_vector: Vec<[ThermodynamicTemperature; NUMBER_OF_NODES]>


}


impl<const N: usize> FluidArray<N> {


    /// obtains a clone of the temperature vector within the CV 
    /// thus obtaining the temperature profile
    pub fn get_temperature_vector(&self) -> Result<
    Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        let mut temperature_vec: Vec<ThermodynamicTemperature> = vec![];

        for temperature in self.temperature_array_current_timestep.iter() {
            temperature_vec.push(*temperature);
        }

        return Ok(temperature_vec);
    }

    /// obtains a clone of the temperature array in Array1 ndarray 
    /// form 
    pub fn get_temperature_array(&self) -> Result< 
    Array1<ThermodynamicTemperature>, ThermalHydraulicsLibError> {

        // converts the fixed sized temperature array (at compile time) 
        // into a dynamically sized ndarray type so we can use solve
        // methods
        let mut temperature_arr: Array1<ThermodynamicTemperature> = 
        Array1::default(N);

        for (idx,temperature) in 
            self.temperature_array_current_timestep.iter().enumerate() {
                temperature_arr[idx] = *temperature;
        }

        Ok(temperature_arr)


    }

    /// obtains a clone of temperature vector, but in reverse format 
    /// this is useful for counter flow heat exchangers 
    pub fn get_reverse_temperature_vector(&self) -> 
    Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        let vec = self.get_temperature_vector()?;

        let reversed_vec = vec.iter().copied().rev().collect();

        Ok(reversed_vec)
    }

    /// sets the temperature vector to a 
    pub fn set_temperature_vector(&mut self,
    temperature_vec: Vec<ThermodynamicTemperature>) -> Result<(), ThermalHydraulicsLibError>{

        let number_of_temperature_nodes = N;

        // check if temperature_vec has the correct number_of_temperature_nodes

        if temperature_vec.len() !=  number_of_temperature_nodes {
            let shape_error = ShapeError::from_kind(
                ErrorKind::IncompatibleShape
            );

            let linalg_error = LinalgError::Shape(shape_error);

            return Err(ThermalHydraulicsLibError::LinalgError
                (linalg_error));

        }

        for (index,temperature) in 
            self.temperature_array_current_timestep.iter_mut().enumerate() {
            *temperature = temperature_vec[index];
        }

        // we also need to ensure that the front and end nodes are 
        // properly synchronised in terms of temperature
        //

        let back_cv_temperature: ThermodynamicTemperature 
        = temperature_vec[0];

        let front_cv_temperature: ThermodynamicTemperature 
        = *temperature_vec.last().unwrap();


        // update enthalpies of control volumes withing

        let material = self.material_control_volume;
        let pressure = self.pressure_control_volume;

        let back_cv_enthalpy = specific_enthalpy(
            material,
            back_cv_temperature,
            pressure
        )?;

        let front_cv_enthalpy = specific_enthalpy(
            material,
            front_cv_temperature,
            pressure
        )?;

        self.back_single_cv.current_timestep_control_volume_specific_enthalpy
            = back_cv_enthalpy;

        self.front_single_cv.current_timestep_control_volume_specific_enthalpy
            = front_cv_enthalpy;


        
        Ok(())
    }

    /// obtains a clone of the temperature array in Array1 ndarray 
    /// form 
    pub fn set_temperature_array(&mut self,
    temperature_arr: Array1<ThermodynamicTemperature>) -> Result<(),
    ThermalHydraulicsLibError> {

        // we'll convert the temperature array into vector form 
        // and use the existing method 
        let mut temperature_vec: Vec<ThermodynamicTemperature> 
        = vec![];

        for temperature in temperature_arr.iter() {
            temperature_vec.push(*temperature)
        }

        self.set_temperature_vector(temperature_vec)


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


/// Contains functions specifically pertaining to fluid trait 
/// implementation 
pub mod fluid_component_trait;
pub use fluid_component_trait::*;

