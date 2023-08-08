use super::lumped_nuclear_structure_inspired_functions::*;
use std::f64::consts::PI;

use ndarray::*;
use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::power::watt;
use uom::si::pressure::atmosphere;
use uom::si::time::second;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::calculations::UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS;
use crate::heat_transfer_lib::
thermophysical_properties::specific_heat_capacity::specific_heat_capacity;
use crate::heat_transfer_lib::
thermophysical_properties::thermal_diffusivity::thermal_diffusivity;
use crate::heat_transfer_lib::
thermophysical_properties::density::density;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::Material;

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
#[derive(Debug,Clone,PartialEq)]
pub struct CartesianConduction1DArray {
    inner_single_cv: SingleCVNode,
    outer_single_cv: SingleCVNode,

    // number of ADDITIONAL nodes within the array 
    // in addition to the inner and outer single cvs
    inner_nodes: usize,

    total_volume: Volume,

    // conductance array
    conductance_array: Array1<ThermalConductance>,

    // volume fraction array 
    volume_fraction_array: Array1<f64>,

    // temperature array current timestep 
    temperature_array_current_timestep: Array1<ThermodynamicTemperature>,

    // temperature_array_next timestep 
    temperature_array_next_timestep: Array1<ThermodynamicTemperature>,

    // VolumetricHeatCapacity array 
    volumetric_heat_capacity_array: Array1<VolumetricHeatCapacity>,
}

impl CartesianConduction1DArray {

    /// calculates the temperature array for the next timestep 
    /// and updates the temperatures and enthalpies of current timestep 
    /// to be that of the next timestep
    pub fn advance_timestep(
        &mut self, timestep: Time) -> Result<(), String>{
        
        let inner_single_cv_mut_reference = &mut self.inner_single_cv;
        let outer_single_cv_mut_reference = &mut self.outer_single_cv;
        let inner_nodes = self.inner_nodes;
        let total_volume = self.total_volume;

        // this is the same as nodesNumber in the  
        // genfoam code, but I'm reformatting it to make it snake case
        let total_number_of_nodes = inner_nodes + 2;

        // heat supplied to the array is always zero
        let q: Power = Power::new::<watt>(0.0);
        // and the power fraction is always zero too 
        // I used snake case here to rename qFraction as q_fraction
        let q_fraction:Array1<f64> = Array::zeros(total_number_of_nodes);

        // probably want to use a function to obtain the conductances 
        // here
        let conductance_array_mut_reference = &mut self.conductance_array;
        let volume_fraction_array_reference = &self.volume_fraction_array;

        // temperature array 
        let temperature_array_current_timestep_reference = 
        &self.temperature_array_current_timestep;

        // rhoCp array 
        let volumetric_heat_capacity_array_reference = 
        &self.volumetric_heat_capacity_array.clone();

        let new_temperature_array: Array1<ThermodynamicTemperature> = 
        advance_timestep_for_specified_conductance_array_cv(
            inner_single_cv_mut_reference,
            outer_single_cv_mut_reference,
            inner_nodes,
            timestep,
            total_volume,
            q,
            &q_fraction,
            temperature_array_current_timestep_reference,
            conductance_array_mut_reference,
            volume_fraction_array_reference,
            volumetric_heat_capacity_array_reference
        ).unwrap();

        self.temperature_array_next_timestep = 
            new_temperature_array.clone();

        // Todo: probably need to synchronise error types in future
        //

        // I'm calculating the inner and outer CV's new enthalpy
        // at the current timestep 
        //
        // This code block deals with setting temperature and 
        // enthalpies for the two nested control volumes
        {
            let inner_node_enthalpy_next_timestep: AvailableEnergy = 
            specific_enthalpy(
                self.inner_single_cv.material_control_volume,
                new_temperature_array[0],
                self.inner_single_cv.pressure_control_volume).unwrap();

            self.inner_single_cv.current_timestep_control_volume_specific_enthalpy 
                = inner_node_enthalpy_next_timestep;

            // do the same for the outer node
            let outer_node_enthalpy_next_timestep: AvailableEnergy = 
            specific_enthalpy(
                self.outer_single_cv.material_control_volume,
                new_temperature_array[total_number_of_nodes-1],
                self.outer_single_cv.pressure_control_volume).unwrap();

            self.outer_single_cv.current_timestep_control_volume_specific_enthalpy 
                = outer_node_enthalpy_next_timestep;

            // I also need to update the TOld vector 
            // This will ensure that the current temperature of the single 
            // cv node is equal to that of the matrix

            // technically speaking, this isn't necessary here, as I'm already 
            // returning the new temperature vector
            // I'll set it explicitly in the advance_timestep function
            //
            //*TOld = T.mapv(
            //    |temperature_value| {
            //        return temperature_value;
            //    }
            //);
            // set liquid cv mass for both cvs,
            // and clear out both the rate_enthalpy_change_vector  
            // and the max_timestep_vector
            // probably also need to update error types in future
            // so I don't keep using unwrap
            self.inner_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
            self.inner_single_cv.rate_enthalpy_change_vector.clear();
            self.inner_single_cv.max_timestep_vector.clear();

            self.outer_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
            self.outer_single_cv.rate_enthalpy_change_vector.clear();
            self.outer_single_cv.max_timestep_vector.clear();
        }

        self.temperature_array_current_timestep = new_temperature_array;
        Ok(())
    }
}
