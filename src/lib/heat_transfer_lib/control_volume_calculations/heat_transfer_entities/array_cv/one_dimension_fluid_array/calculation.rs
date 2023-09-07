use ndarray::*;
use ndarray_linalg::error::LinalgError;
use uom::num_traits::Zero;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;

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

use super::FluidArray;

/// this implementation deals with lateral connections 
///
/// the convention is to supply an average conductance 
/// as well as a temperature array
///
/// at the end of the connection phase, one can then use 
/// the advance_timestep method to calculate the new 
/// temperature array
impl<const NUMBER_OF_NODES: usize> FluidArray<NUMBER_OF_NODES>{
    
    /// advances timestep for the fluid array 
    /// given a fixed timestep
    pub fn advance_timestep_with_mass_flowrate(&mut self,
        timestep: Time,
        mass_flowrate: MassRate) 
        -> Result<(), ThermalHydraulicsLibError>{

        self.mass_flowrate = mass_flowrate;
        //advance_timestep_fluid_node_array_pipe_high_peclet_number(
        //    back_cv_ptr_in_loop.deref_mut(),
        //    front_cv_ptr_in_loop.deref_mut(),
        //    number_of_nodes,
        //    timestep,
        //    total_volume,
        //    heater_steady_state_power,
        //    steel_temp_at_present_timestep_ptr_in_loop.deref_mut(),
        //    &mut conductance_vector,
        //    fluid_temp_vec_ptr_in_loop.deref_mut(),
        //    therminol_mass_flowrate,
        //    fluid_vol_fraction_ptr_in_loop.deref_mut(),
        //    &mut fluid_rho_cp_array,
        //    q_fraction_ptr_in_loop.deref_mut(),
        //)?;

        // there will always be at least 2 nodes

        if NUMBER_OF_NODES <= 1 {
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
        // The old fluid temperature will need to be used to calculate 
        // new specific enthalpy for the system

        let mut temperature_vector: Array1<ThermodynamicTemperature> = 
        self.get_temperature_array()?.map(
            |temp_kelvin_ptr: &ThermodynamicTemperature| {

                return *temp_kelvin_ptr;
            }

        );
        // now let's start calculation 
        //
        if NUMBER_OF_NODES > 1 {
            let mut coefficient_matrix: Array2<ThermalConductance> = 
            Array::zeros((NUMBER_OF_NODES, NUMBER_OF_NODES));

            let mut power_source_vector: 
            Array1<Power> = Array::zeros(NUMBER_OF_NODES);

            // ascertain if we have forward flow 

            let forward_flow: bool = self.
                mass_flowrate.ge(&MassRate::zero());

            // back node calculation (first node)
            {
                let volume_fraction_array: Array1<f64> = 
                self.volume_fraction_array.iter().map(
                    |&vol_frac| {
                        vol_frac
                    }
                ).collect();

                // rho_cp is determined by bulk temperature 
                //

                // for the first node, also called the back node
                // energy balance is: 
                // m c_p dT/dt = -H (T - T_solid) - m_flow h_fluid(T) 
                // + m_flow h_fluid(adjacent T) + q
                // of all these terms, only the m cp dT/dt term and HT 
                // 
                // We separate this out to get:
                //
                // m cp T / dt + HT = 
                //
                // HT_solid - m_flow h_fluid(T_old) 
                // + m_flow h_fluid(adjacent T_old) + m cp / dt (Told)
                // + q
                //
                // belong in the M matrix, the rest belong in S
                //coefficient_matrix[[0,0]] = volume_fraction_array[0] * rho_cp[0] 
                //    * total_volume / dt + solid_fluid_conductance_array[0];

                //// the first part of the source term deals with 
                //// the flow direction independent terms

                //let h_fluid_last_timestep: AvailableEnergy = 
                //back_single_cv.current_timestep_control_volume_specific_enthalpy;

                //// now this makes the scheme semi implicit, and we should then 
                //// treat the scheme as explicit

                //power_source_vector[0] = solid_fluid_conductance_array[0] *
                //    last_timestep_temperature_solid[0] 
                //    - mass_flowrate * h_fluid_last_timestep 
                //    + last_timestep_temperature_fluid[0] * total_volume * 
                //    volume_fraction_array[0] * rho_cp[0] / dt 
                //    + q * q_fraction[0] 
                //    + total_enthalpy_rate_change_back_node ;

                //// the next part deals with the inflow
                //// m_flow h_fluid(adjacent T_old)
                ////
                //// now if the advection interaction is done correctly, 
                ////
                //// (advection) ----- (back cv) --------> fwd
                ////
                //// then in a frontal flow condition, the enthalpy flows in 
                //// would already have been accounted for
                ////
                //// but in the case of backflow, then fluid from the node
                //// in front will flow into this fluid node 
                //// that is node 1 

                //// so if mass flowrate is <= 0 , then we will calculate 
                //// backflow conditions

                //if !forward_flow {
                //    // first, get enthalpy of the node in front 

                //    let enthalpy_of_adjacent_node_to_the_front: AvailableEnergy = 
                //    specific_enthalpy(
                //        back_single_cv.material_control_volume,
                //        last_timestep_temperature_fluid[1],
                //        back_single_cv.pressure_control_volume).unwrap();

                //    // now if mass flowrate is less than zero, then 
                //    // we receive enthalpy from the front cv 
                //    //
                //    // But we need to subtract the negative mass flow if 
                //    // that makes sense, or at least make it absolute
                //    // - (-m) * h = m * h
                //    //

                //    power_source_vector[0] += 
                //    mass_flowrate.abs() * enthalpy_of_adjacent_node_to_the_front;

                //    // additionally, in backflow situations, the mass 
                //    // flow out of this cv is already accounted for 
                //    // so don't double count 

                //    power_source_vector[0] += 
                //    mass_flowrate * h_fluid_last_timestep;



                //}
            }
        }

        self.clear_vectors()?;
        todo!()
    }

    /// clears all vectors for next timestep
    /// This is important for the advance timestep method
    pub fn clear_vectors(&mut self) 
    -> Result<(), ThermalHydraulicsLibError>{

        self.lateral_adjacent_array_conductance_vector.clear();
        self.lateral_adjacent_array_temperature_vector.clear();
        Ok(())
    }
}
