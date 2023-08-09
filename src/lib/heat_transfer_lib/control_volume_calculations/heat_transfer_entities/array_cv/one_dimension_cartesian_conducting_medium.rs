use super::lumped_nuclear_structure_inspired_functions::*;
use std::f64::consts::PI;

use ndarray::*;
use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::power::watt;
use uom::si::pressure::atmosphere;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::time::second;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::single_cv_node;
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


    /// constructs a new instance of the CartesianConduction1DArray
    pub fn new(material: Material,
    initial_uniform_temperature: ThermodynamicTemperature,
    uniform_pressure: Pressure,
    inner_nodes: usize,
    total_length: Length) -> Result<HeatTransferEntity, String> {
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
            uniform_pressure)?;

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

        return Ok(heat_transfer_entity);
    }

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
    /// The temperatures of T[0] and T[n-1] correspond to surface 
    /// temperatures
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

    /// constructor method which helps set the volume fraction array
    /// This returns a volume fraction array
    /// based on the following diagram:
    ///
    ///
    ///
    /// Tmax            T[0]           T[1]          T[n-1]         T_ambient
    /// [ignored]       T_innersingleCV             T_outersingleCV [ignored]
    ///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
    ///        [ignored]                                    [ignore]
    /// 
    ///
    /// Between each temperature node, there is a thermal resistor 
    /// each resistor has resistance equivalent to delta x 
    /// long 
    ///
    /// where delta x = total length / number of resistors

    /// the volume fraction rray helps to determine the thermal inertia 
    /// of each node 
    /// now each node would have a thermal inertia corresponding to 
    /// a thickness of delta_x 
    ///
    /// However, the nodes at both boundaries would have a length 
    /// equivalent to half of delta_x because they are at the boundaries
    /// These would constitute half the thermal inertia of a typical node 
    fn construct_volume_fraction_array(&mut self) 
    -> Result<Array1<f64>, String> {

        todo!("recompute thermal inertia for the volumes");

        let number_of_temperature_nodes: usize = 
        self.inner_nodes + 2;

        let number_of_thermal_resistors: usize = 
        number_of_temperature_nodes - 1;

        let delta_x: Length = self.total_length/
        number_of_thermal_resistors as f64;

        // now we can compute the volume fraction based on length 
        // because the basis cross sectional area is the same 
        // throughout

        let volume_fraction_per_inner_node: f64 = 
        (delta_x/self.total_length).value;

        let volume_fraction_boundary_node: f64 
        = volume_fraction_per_inner_node * 0.5;

        let mut volume_fraction_array: Array1<f64> = 
        Array::zeros(number_of_temperature_nodes);

        // set all volume fractions to the inner node volume fractions
        volume_fraction_array.fill(volume_fraction_per_inner_node);

        // set the boundary nodes to half that fraction 
        volume_fraction_array[0] = volume_fraction_boundary_node;
        volume_fraction_array[number_of_temperature_nodes-1] = 
            volume_fraction_boundary_node;

        // assert if they add up to 1.0

        let vol_fraction_sum: f64 = volume_fraction_array.sum();
        assert_eq!(1.0, vol_fraction_sum);


        return Ok(volume_fraction_array);
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

    #[inline]
    pub fn get_bulk_temperature(&mut self) -> 
    Result<ThermodynamicTemperature,String>{

        // for now, doing it quick and dirty, i'm going to obtain a volume 
        // averaged temperature 

        let volume_fraction_array_reference = 
        &self.volume_fraction_array;
        let temperature_array_reference = 
        &self.temperature_array_current_timestep;

        let mut vol_averaged_temperature_array_values: Array1<f64> 
        = Array::default(temperature_array_reference.len());

        for (idx, temperature_reference) in 
            temperature_array_reference.iter().enumerate() {
                //get the vol fraction 

                let vol_fraction: f64 = 
                volume_fraction_array_reference[idx];

                let vol_avg_temperature_component: f64
                = vol_fraction * (temperature_reference.get::<kelvin>());

                vol_averaged_temperature_array_values[idx] = 
                    vol_avg_temperature_component;

            }

        // sum it all up (these are float values) 

        let vol_averaged_temperature_kelvin: f64 
        = vol_averaged_temperature_array_values.sum();

        return Ok(ThermodynamicTemperature::new
            ::<kelvin>(vol_averaged_temperature_kelvin));


    }
}
