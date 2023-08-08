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
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
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
use crate::heat_transfer_lib::thermophysical_properties::volumetric_heat_capacity::rho_cp;

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
    
    // total length for the 1D array
    total_length: Length,

    // volume fraction array 
    volume_fraction_array: Array1<f64>,

    // temperature array current timestep 
    temperature_array_current_timestep: Array1<ThermodynamicTemperature>,

    // temperature_array_next timestep 
    temperature_array_next_timestep: Array1<ThermodynamicTemperature>,

    /// control volume material 
    material_control_volume: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,
}

impl CartesianConduction1DArray {

    fn get_current_timestep_rho_cp(
        material: Material,
        temperature_array_current_timestep_reference: &Array1<ThermodynamicTemperature>,
        pressure: Pressure,
    ) -> 
    Result<Array1<VolumetricHeatCapacity>,String>{

        let rho_cp_array:Array1<VolumetricHeatCapacity> 
        = temperature_array_current_timestep_reference.map(
            // here is the closure (function) i shall use
            |temperature_reference: &ThermodynamicTemperature|{
                let temperature = *temperature_reference;

                let rho_cp: VolumetricHeatCapacity 
                = rho_cp( material, temperature, pressure).unwrap();

                return rho_cp;
            }
        );

        return Ok(rho_cp_array);

    }

    /// This returns a conductance array
    /// based on the following diagram:
    ///
    ///
    ///
    /// Tmax            T[0]           T[1]          T[n-1]         T_ambient
    /// [ignored]       T_innersingleCV             T_outersingleCV [ignored]
    ///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
    ///        [ignored]                                    [ignore]
    ///
    /// Basically, H[0] and H[n] are first and last elements of the 
    /// conductance arrays, these values are ignored because I don't 
    /// want to rewrite code for the advance_timestep function as  
    /// far as possible
    ///
    /// to construct the array, we consider that we have n nodes
    /// connected via thermal resistors 
    ///
    /// for thermal conductivity in 1D cartesian coordinates ,
    /// the thermal conductance is kA/L 
    ///
    /// Where A is the cross sectional area, and L is the distance 
    /// from node to node
    ///
    /// A is determined for 1D cross sectional area using a global constant
    ///
    /// For our case, we assume that the conductance is contributed to 
    /// by two thermal resistances.
    ///
    /// For example H[1] is contributed to by kA/(0.5L) for T[0] and
    /// kA/(0.5L) for T[1]
    ///
    ///
    ///
    fn get_current_timestep_conductance_array(
        material: Material,
        temperature_array_curnent_timestep_reference: &Array1<ThermodynamicTemperature>,
        pressure: Pressure,
        total_length: Length,
    ) -> 
    Result<Array1<ThermalConductance>,String> {

        // for 1D calcs, we use a basis area
        let basis_area: Area = Area::new::<square_meter>( 
            UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS);

        // we can obtain a map of thermal conductivity first 

        let thermal_conductivity_array: Array1<ThermalConductivity> 
        = temperature_array_curnent_timestep_reference.map(
            |temperature_reference: &ThermodynamicTemperature| {

                let temperature = *temperature_reference;

                let k: ThermalConductivity 
                = thermal_conductivity( 
                    material, temperature, pressure).unwrap();

                return k;
            }
        );
        // now we are dealing with a number of thermal resistors 

        let number_of_thermal_resistors: usize = 
        thermal_conductivity_array.len() - 1;

        let delta_x: Length = total_length/
        number_of_thermal_resistors as f64;

        let mut thermal_conductance_vector: Vec<ThermalConductance> 
        = vec![];

        for index in 0..number_of_thermal_resistors {

            // get first thermal resistance at the index
            let thermal_conductance_idx: ThermalConductance
            = thermal_conductivity_array[index] * basis_area/
            delta_x;

            let thermal_conductance_idx_plus_one: ThermalConductance 
            = thermal_conductivity_array[index + 1] * basis_area/
            delta_x;

            // sum both thermal resistance 

            let total_thermal_resistance = 
            1.0/thermal_conductance_idx + 
            1.0/thermal_conductance_idx_plus_one;

            // final thermal conductance 

            let thermal_conductance_idx_and_idx_plus_one = 
            1.0/total_thermal_resistance;

            // push to conductance vector
            thermal_conductance_vector.push(
            thermal_conductance_idx_and_idx_plus_one);
        }
        
        // now let's make the array with number of thermal 
        // resistors plus 2, the first and last elements 
        // are zero because I want to save effort rewriting the 
        // finite difference code

        let mut thermal_conductance_array: Array1<ThermalConductance> 
        = Array::zeros(number_of_thermal_resistors + 2);

        for (idx, conductance_reference) in 
            thermal_conductance_vector.iter().enumerate() {

                // load conductance array at index plus 1 

                thermal_conductance_array[idx + 1] 
                = *conductance_reference;

        }
        // ok done!
        Ok(thermal_conductance_array)
    }

    /// calculates the temperature array for the next timestep 
    /// and updates the temperatures and enthalpies of current timestep 
    /// to be that of the next timestep
    pub fn advance_timestep(
        &mut self, timestep: Time) -> Result<(), String>{
        
        let inner_single_cv_mut_reference = &mut self.inner_single_cv;
        let outer_single_cv_mut_reference = &mut self.outer_single_cv;
        let inner_nodes = self.inner_nodes;

        // for 1D calcs, we use a basis area
        let basis_area: Area = Area::new::<square_meter>( 
            UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS);
        let total_volume = self.total_length*basis_area;

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
        // 
        let material = self.material_control_volume;
        let pressure = self.pressure_control_volume;
        // temperature array 
        let temperature_array_current_timestep_reference = 
        &self.temperature_array_current_timestep;


        let mut conductance_array: Array1<ThermalConductance> = 
        Self::get_current_timestep_conductance_array(
            material,
            temperature_array_current_timestep_reference,
            pressure,
            self.total_length
        )?;
        let conductance_array_mut_reference = &mut conductance_array;
        let volume_fraction_array_reference = &self.volume_fraction_array;


        // rhoCp array 
        // probably get rid of unwrap later

        let volumetric_heat_capacity_array = 
        Self::get_current_timestep_rho_cp(
            material,
            temperature_array_current_timestep_reference,
            pressure
        ).unwrap();

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
            &volumetric_heat_capacity_array
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
