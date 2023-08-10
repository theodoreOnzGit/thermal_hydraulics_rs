use super::lumped_nuclear_structure_inspired_functions::*;
use std::f64::consts::PI;

use approx::assert_relative_eq;
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
    material_control_volume: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,
}

impl CartesianConduction1DArray {

    /// returns a clone of the temperature_array_current_timestep
    pub fn get_temperature_vector(&mut self) -> 
    Result<Vec<ThermodynamicTemperature>, String> {
        let temp_array = self.temperature_array_current_timestep.clone();

        let mut temp_vector: Vec<ThermodynamicTemperature> = vec![];

        for temp_reference in temp_array.iter() {
            temp_vector.push(*temp_reference)
        }

        Ok(temp_vector)
    }
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
    /// The trick is in calculating L, the node to node distance 
    /// for purposes of calculating thermal resistance 
    ///
    /// consider a four node system:
    ///
    /// surf1 [0]    [1]     [2]     [3]    surf2
    /// | --- * ------- x ------ x ------ * --- |
    ///
    /// |<----------- total_length ------------>|
    ///
    /// Now L would be the full length, however 
    /// the total length of the thermal resistors 
    /// would be slightly less than total_length because we are not 
    /// measuring the full thermal resistance from surface 1 to 
    /// surface 2 
    ///
    /// the full length of the thermal resistors would only account 
    /// for node [0] to node [3]
    ///
    /// so we can find the node to node distance assuming the nodes 
    /// are spaced equally apart and the distance from the boundary 
    /// nodes to surfaces are about 0.5L
    ///
    /// surf1 [0]    [1]     [2]     [3]    surf2
    /// | --- * ------- x ------ x ------ * --- |
    ///  0.5L     L         L        L      0.5L
    /// |<----------- total_length ------------>|
    ///
    /// So the total length is 4L, and 4 is the number of temperature 
    /// nodes
    ///
    /// however, each thermal resistor is connected by a resistance of L
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

        let number_of_temperature_nodes: usize = 
        thermal_conductivity_array.len();

        let number_of_thermal_resistors: usize = 
        number_of_temperature_nodes - 1;

        let delta_x: Length = total_length/
        number_of_temperature_nodes as f64;

        let mut thermal_conductance_vector: Vec<ThermalConductance> 
        = vec![];

        for index in 0..number_of_thermal_resistors {

            // get first thermal resistance at the index
            let thermal_conductance_idx: ThermalConductance
            = thermal_conductivity_array[index] * basis_area/
            (0.5*delta_x);

            let thermal_conductance_idx_plus_one: ThermalConductance 
            = thermal_conductivity_array[index + 1] * basis_area/
            (0.5*delta_x);

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
    ///
    /// However, the nodes at both boundaries would have a length 
    /// equivalent to half of delta_x because they are at the boundaries
    /// This would only apply to the thermal resistor
    ///
    /// In terms of thermal inertia, this would not make a difference 
    /// because the last node would have the same thermal inertia as 
    /// any other node. Therefore, length is evenly divided among all 
    /// nodes
    ///
    /// Consider again a four node system with the thermal resistors 
    /// of length L within the bulk and 0.5L at the boundary
    ///
    /// surf1 [0]    [1]     [2]     [3]    surf2
    /// | --- * ------- x ------ x ------ * --- |
    ///  0.5L     L         L        L      0.5L
    /// |<----------- total_length ------------>|
    ///
    /// We must account for the full thermal inertia of the system,
    /// and therefore, the whole length is taken into account. 
    /// The thermal inertia of each array is represented by length L
    ///
    ///
    fn construct_volume_fraction_array(&mut self) 
    -> Result<Array1<f64>, String> {

        let number_of_temperature_nodes: usize = 
        self.inner_nodes + 2;

        let delta_x: Length = self.total_length/
        number_of_temperature_nodes as f64;

        // now we can compute the volume fraction based on length 
        // because the basis cross sectional area is the same 
        // throughout

        let volume_fraction_per_inner_node: f64 = 
        (delta_x/self.total_length).value;

        let mut volume_fraction_array: Array1<f64> = 
        Array::zeros(number_of_temperature_nodes);

        // set all volume fractions to the inner node volume fractions
        volume_fraction_array.fill(volume_fraction_per_inner_node);

        // assert if they add up to 1.0

        let vol_fraction_sum: f64 = volume_fraction_array.sum();
        assert_relative_eq!(
            1.0,
            vol_fraction_sum,
            epsilon = 0.001);

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

    /// gets bulk temperature of the array cv based on volume fraction 
    /// now, for solid and liquid, that would be sort of 
    /// a good approximation since the boussinesq approximation
    /// may work well for liquids
    ///
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

    /// gets the maximum timestep from the one dimensional 
    /// control volume for cartesian conduction
    pub fn get_max_timestep(&mut self, 
    max_temperature_change: TemperatureInterval) -> Result<Time,String>{
        
        // start with an empty timestep vector
        let mut max_timestep_vector: Vec<Time> = vec![];

        // let's find alpha, 
        let control_vol_pressure: Pressure = self.pressure_control_volume;
        let control_vol_material: Material = self.material_control_volume;
        let control_vol_temperature_array: Array1<ThermodynamicTemperature> 
        = self.temperature_array_current_timestep.clone();

        // now let's map the alpha 
        // this is quick and dirty cos i used unwrap

        let thermal_diffusivity_array: Array1<DiffusionCoefficient>
        = control_vol_temperature_array.map(
            |temperature_reference| {

                thermal_diffusivity(control_vol_material, 
                    *temperature_reference, 
                    control_vol_pressure).unwrap()
            }

        );

        // let's calculate the internal lengthscale
        let number_of_temperature_nodes: usize = 
        self.inner_nodes + 2;

        let delta_x: Length = self.total_length/
        number_of_temperature_nodes as f64;

        // now for the array CV, implicit schemes are used. Therefore,
        // the threshold for stability is higher, at 1.0 
        // out of some caution, let me use 0.8
        let max_mesh_fourier_number: f64 = 0.8;

        // for the minimum conduction timescale, we need the 
        // maximum alpha
        //
        // I'm using this closure to find the maximum, rather than 
        // a manual for loop
        let thermal_diffusivity_coeff_opt 
        = thermal_diffusivity_array.iter().max_by(
            |alpha_1, alpha_2| {
                // a and b represent two typical items in the array 
                let a = &alpha_1.value;
                let b = &alpha_2.value;
                a.total_cmp(b)
            });

        let thermal_diffusivity_coeff: DiffusionCoefficient = 
        match thermal_diffusivity_coeff_opt {
            Some(alpha_reference) => *alpha_reference,
            None => {
                // the none case should NOT happen at all, I'm just 
                // otherwise it means that it's impossible to get thermal 
                // diffusivity
                // providing a fallback mechanism
                let bulk_temp =self.get_bulk_temperature().unwrap();

                thermal_diffusivity(control_vol_material, 
                    bulk_temp, 
                    control_vol_pressure).unwrap()
            },
        };

        // timescales for conduction of this array
        let max_conduction_timescale: Time = max_mesh_fourier_number * 
        delta_x *
        delta_x / 
        thermal_diffusivity_coeff;

        max_timestep_vector.push(max_conduction_timescale);

        // we also need to take into account the timescales of each of 
        // the control volume at the boundaries

        let max_timestep_front_cv: Time = 
        self.outer_single_cv.get_max_timestep(max_temperature_change)?;

        let max_timestep_back_cv: Time = 
        self.inner_single_cv.get_max_timestep(max_temperature_change)?;

        max_timestep_vector.push(max_timestep_front_cv);
        max_timestep_vector.push(max_timestep_back_cv);

        // now we just find the minimum in this max_timestep_vector 
        //
        // this would help to find the timestep vector without using 
        // a manual for loop, though of course, the for loop is more 
        // readable for rust beginners

        //let mut maximum_timestep = Time::new::<second>(75_f64);

        //for time in max_timestep_vector.iter() {
        //    if *time < maximum_timestep {
        //        maximum_timestep = *time;
        //    }
        //} 
        //
        // will probably want to watch the edge BC timesteps as well

        let maximum_timestep: Time = 
        *max_timestep_vector.iter().min_by(
            |time_1, time_2| {
                // a and b represent two typical items in the array 
                let a = &time_1.value;
                let b = &time_2.value;
                a.partial_cmp(b).unwrap()
            }).unwrap();

        // all right done!

        return Ok(maximum_timestep);
    }
}
