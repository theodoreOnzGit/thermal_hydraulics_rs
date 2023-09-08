use super::SolidColumn;
use uom::{si::{f64::*, thermodynamic_temperature::kelvin, temperature_interval::degree_celsius}, num_traits::Zero};
use crate::{thermal_hydraulics_error::ThermalHydraulicsLibError, heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity};
use ndarray_linalg::error::LinalgError;
use ndarray::*;
use uom::si::power::watt;
use crate::heat_transfer_lib::thermophysical_properties::volumetric_heat_capacity::rho_cp;

/// This deals with the calculations of the solid column array
/// at the end of the connection phase, one can then use 
/// the advance_timestep method to calculate the new 
/// temperature array
impl SolidColumn {

    /// advances timestep for the solid array column
    /// given a fixed timestep
    pub fn advance_timestep(&mut self,
        timestep: Time,) 
        -> Result<(), ThermalHydraulicsLibError>{

        // there must always be at least 2 nodes

        let number_of_nodes = self.len();

        if number_of_nodes <= 1 {
            return Err(LinalgError::Shape(
                ShapeError::from_kind(
                    ErrorKind::OutOfBounds
                )).into());
        }
        // First things first, we need to set up 
        // how the CV interacts with the internal array
        // here is heat added to CV

        let back_cv_rate_enthalpy_change_vector: Vec<Power> = 
        self.back_single_cv.rate_enthalpy_change_vector.clone();

        let front_cv_rate_enthalpy_change_vector: Vec<Power> = 
        self.front_single_cv.rate_enthalpy_change_vector.clone();

        // compute power source for back node

        let mut total_enthalpy_rate_change_back_node = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            back_cv_rate_enthalpy_change_vector.iter() {

                total_enthalpy_rate_change_back_node += *enthalpy_chg_rate;
            }

        // then the front node,

        let mut total_enthalpy_rate_change_front_node = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            front_cv_rate_enthalpy_change_vector.iter() {

                total_enthalpy_rate_change_front_node += *enthalpy_chg_rate;
            }


        // this front and back nodes will be an extra term added to the 
        // heat source vector S
        //
        // The old solid temperature will need to be used to calculate 
        // new specific enthalpy for the system
        //
        // We need a blank temperature array to start the iteration 
        // process, so we just copy the old temperature array over 
        // as the initial guess (i think?)

        let mut new_temperature_array: Array1<ThermodynamicTemperature> = 
        self.get_temperature_array()?.map(
            |temp_kelvin_ptr: &ThermodynamicTemperature| {

                return *temp_kelvin_ptr;
            }

        );
        // now let's start calculation 
        //
        let mut coefficient_matrix: Array2<ThermalConductance> = 
        Array::zeros((number_of_nodes, number_of_nodes));

        let mut power_source_vector: 
        Array1<Power> = Array::zeros(number_of_nodes);

        // obtain some important parameters for calculation
        let material = self.material_control_volume;
        let pressure = self.pressure_control_volume;
        let bulk_temperature = self.get_bulk_temperature()?;
        let total_volume = self.total_length *  self.xs_area;
        let dt = timestep;
        let node_length = self.total_length / number_of_nodes as f64;
        // for the fluid array, we do not keep track of node enthalpies,
        // instead, we keep track of temperatures and then calculate 
        // enthalpy changes using rho_cp calculated at the current 
        // timestep
        let volume_fraction_array: Array1<f64> = 
        self.volume_fraction_array.iter().map(
            |&vol_frac| {
                vol_frac
            }
        ).collect();
        // rho_cp is determined by the temperature array 
        let rho_cp: Array1<VolumetricHeatCapacity> = 
        self.temperature_array_current_timestep.iter().map(
            |&temperature| {
                rho_cp(material, temperature, pressure).unwrap()
            }
        ).collect();


        // need to determine the sort of connections 
        // we assume axial conduction always occurs unless it can be 
        // neglected
        let lateral_power_sources_connected: bool 
        = self.q_vector.len() > 0; 
        let lateral_temperature_arary_connected: bool 
        = self.lateral_adjacent_array_conductance_vector.len() > 0;

        let axial_conduction_only: bool = !lateral_power_sources_connected
        && !lateral_temperature_arary_connected;

        // if there is lateral conduction, then construct matrices 
        // which take that into account 
        let mut sum_of_lateral_conductances: Array1<ThermalConductance>
        = Array1::zeros(number_of_nodes);

        // conductances will need to be summed over each node 
        //
        // i will also need to make sure that there are actually 
        // lateral connections to this array in the first place 
        // though!


        if lateral_temperature_arary_connected {

            for (idx, node_conductance) in 
                sum_of_lateral_conductances.iter_mut().enumerate() {

                    // I will need to index into the array
                    // and sum the conductances 

                    for lateral_conductance in 
                        &self.lateral_adjacent_array_conductance_vector[idx] {
                            *node_conductance += *lateral_conductance;
                        }

                }
        }

        // we need to do the same for the q and q fractions
        let mut sum_of_lateral_power_sources: Array1<Power>
        = Array::zeros(number_of_nodes);

        if lateral_power_sources_connected {
            // again index into each node, multiply q by the 
            // q fraction 
            //
            // suppose there are two power sources in a three node 
            // system 
            //
            // 
            // [P1a, P1b, P1c] 
            // [P2a, P2b, P2c] 
            //
            // --> these are power ndarray
            //
            // The total power contributed to each node is:
            //
            //
            // [Pa, Pb, Pc]
            //
            // --> this is the sum_of_lateral_power_sources
            //
            // The simplest way would be to multiply the power 
            // q by its respective fractions to obtain a vector
            // of power arrays

            let mut power_ndarray_vector: Vec<Array1<Power>>
            = vec![];

            for (power_source_idx, q_reference) in self.q_vector.iter().enumerate() {

                // multiply q by qfraction


                let power_frac_array: Array1<f64>
                = self.q_fraction_vector[power_source_idx].clone();

                let power_ndarray: Array1<Power>
                = power_frac_array.map(
                    |&power_frac| {
                        power_frac * (*q_reference)
                    }

                );


                power_ndarray_vector.push(power_ndarray);

            }

            // this part construts
            // the sum_of_lateral_power_sources
            // vector
            //
            // [Pa, Pb, Pc]
            //
            // at first, it's just a vector of zero power 
            // 
            // [0, 0, 0]
            //
            // The outer loop gets the node index, 
            // for example, node 1, 
            //
            // and obtains a reference to the node power source,
            // 0 watts
            //

            for (node_idx, node_power_source) in 
                sum_of_lateral_power_sources.iter_mut().enumerate() {

                    // at the node index i, and 
                    // with a reference to the node power source,
                    // I want to sum all the power sources
                    //
                    // So I start indexing into each power vector

                    for power_vector in power_ndarray_vector.iter(){

                        // from each power vector, I pull out the 
                        // relevant node power at the correct 
                        // node

                        let node_power_contribution = 
                        power_vector[node_idx];

                        // then I add it to the power source

                        *node_power_source 
                        += node_power_contribution;


                    }
                    // end inner for loop

                }
            // end sum of lateral sources for loop

        }
        // end if for lateral_power_sources_connected


        // now that we've gotten all the important properties, we can 
        // start matrix construction


        // back node calculation (first node)
        {

            // for the first node, also called the back node
            // energy balance is: 
            // m c_p dT/dt = -\sum H (T - T_lateral) + q
            // of all these terms, only the m cp dT/dt term and HT 
            // 
            // We separate this out to get:
            //
            // m cp T / dt + \sum HT = 
            // \sum HT_lateral + m cp / dt (Told)
            // + q
            //
            // so we will need to determine sum H and sum HT_lateral
            // as sum H is the relevant coefficient in the coefficient_matrix




            // Now I'm ready to construct the M matrix
            // belong in the M matrix, the rest belong in S
            coefficient_matrix[[0,0]] = 
                volume_fraction_array[0] * rho_cp[0] 
                * total_volume / dt + sum_of_lateral_conductances[0];


            // now this makes the scheme semi implicit, and we should then 
            // treat the scheme as explicit

            todo!("lateral temperature arrays need to be connected");
            power_source_vector[0] = sum_of_lateral_conductances[0] *
                self.temperature_array_current_timestep[0] 
                + self.temperature_array_current_timestep[0] * total_volume * 
                volume_fraction_array[0] * rho_cp[0] / dt 
                + sum_of_lateral_power_sources[0]
                + total_enthalpy_rate_change_back_node ;



        }
        // end back node calculation code block

        
        // bulk node calculations 
        // only takes into account lateral conductances
        if number_of_nodes > 2 {
            // loop over all nodes from 1 to n-2 (n-1 is not included)
            for i in 1..number_of_nodes-1 {
                // the coefficient matrix is pretty much the same, 
                // we only consider solid fluid conduction, and the 
                // thermal inertia terms

                coefficient_matrix[[i,i]] = volume_fraction_array[i] * rho_cp[i] 
                    * total_volume / dt + sum_of_lateral_conductances[i];

                // we also consider outflow using previous timestep 
                // temperature, 
                // assume back cv and front cv material are the same

            todo!("lateral temperature arrays need to be connected");
                // basically, all the power terms remain 
                power_source_vector[i] = sum_of_lateral_conductances[i] *
                    self.temperature_array_current_timestep[i] 
                    + self.temperature_array_current_timestep[i] * total_volume * 
                    volume_fraction_array[i] * rho_cp[i] / dt 
                    + sum_of_lateral_power_sources[i];


            }
            // end for loop for bulk nodes
        }
        // end bulk node calculation block

        // front node (last node) calculation 
        {
            // we still follow the same guideline
            // m cp T / dt + HT = 
            //
            // HT_solid + m cp / dt (Told)
            // + q
            let i = number_of_nodes-1;

            coefficient_matrix[[i,i]] = volume_fraction_array[i] * rho_cp[i] 
                * total_volume / dt + sum_of_lateral_conductances[i];
            // the first part of the source term deals with 
            // the flow direction independent terms

            // now this makes the scheme semi implicit, and we should then 
            // treat the scheme as explicit
            //
            // now I should subtract the enthalpy outflow from the 
            // power source vector here, 
            // but if done correctly, it should already be accounted 
            // for in total_enthalpy_rate_change_front_node
            //
            // so i shouldn't double count

            todo!("lateral temperature arrays need to be connected");
            power_source_vector[i] = sum_of_lateral_conductances[i] *
                self.temperature_array_current_timestep[i] 
                + self.temperature_array_current_timestep[i] * total_volume * 
                volume_fraction_array[i] * rho_cp[i] / dt 
                + sum_of_lateral_power_sources[i] 
                + total_enthalpy_rate_change_front_node ;


        }

        // the above takes care of lateral conduction and thermal 
        // inertia, if there was no need to worry about axial conduction 
        // then we can just follow through
        //
        // We can neglect axial conduction only if the lateral 
        // conduction is much greater than axial conduction 
        let neglect_axial_conduction: bool = {

            // how shall we know if the axial conduction is to be 
            // neglected?


            let axial_power_scale_insignificant = || -> bool {

                // we can calculate a typical power scale for axial 
                // conduction and compare it to the radial conduction

                // first calculate axial conduction thermal 
                // resistance

                let average_thermal_conductivity = 
                thermal_conductivity(
                    material,
                    bulk_temperature,
                    pressure,
                ).unwrap();

                // now calculate axial thermal conductance
                let average_axial_thermal_conductance: ThermalConductance 
                = average_thermal_conductivity * self.xs_area 
                / node_length;

                // the max temperature gradient is the max temperature 
                // minus the min temperature 

                let temp_value_kelvin_integer_vector: Array1<u32> = 
                self.temperature_array_current_timestep.map(
                    |&temperature|{
                        let temp_value_kelvin = temperature.get::<kelvin>();

                        temp_value_kelvin.ceil() as u32
                    }
                );

                let max_temp_val_kelvin: u32 = 
                *temp_value_kelvin_integer_vector.iter().max().unwrap();

                let min_temp_val_kelvin: u32 = 
                *temp_value_kelvin_integer_vector.iter().min().unwrap();

                let approx_axial_temp_diff_val_kelvin: f64 = 
                max_temp_val_kelvin as f64 
                - min_temp_val_kelvin as f64;

                let axial_power_scale: Power = 
                TemperatureInterval::new::<degree_celsius>(
                    approx_axial_temp_diff_val_kelvin)
                * average_axial_thermal_conductance;
                

                // we need to compare this against the radial temperature 
                // differences
                //
                // as well as power inputs 
                //
                // unfortunately, the temperatures are nested in a 
                // vector of arrays or in essence a 2D array

                let mut lateral_power_sum: Power = Power::zero();

                for &power in &self.q_vector {
                    lateral_power_sum += power;
                }


                // now an estimate for the radial thermal conductance 

                // I'm going to take the sum of lateral conductances 
                // and take that as some average 

                let sum_of_lateral_conductance_estimate: ThermalConductance 
                = sum_of_lateral_conductances[0];
                





                true
            };



            // first let's deal with the case that it's axial conduction 
            // only 

            if axial_conduction_only {
                false
            } else {
                true
            }

        };
        // end code block for neglecting axial conduction


        
        
        todo!()
    }
}
