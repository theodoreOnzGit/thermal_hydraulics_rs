use super::NonInsulatedParallelFluidComponent;
use uom::si::f64::*;
use uom::si::power::watt;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::volumetric_heat_capacity::try_get_rho_cp;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use std::thread::JoinHandle;
use std::thread;
use ndarray::*;
use ndarray_linalg::error::LinalgError;
use uom::num_traits::Zero;

impl NonInsulatedParallelFluidComponent {

    /// advances timestep for each HeatTransferEntity within the 
    /// NonInsulatedPipe
    /// treats the pipe as a single tube
    #[inline]
    fn _advance_timestep_single_tube(&mut self, 
    timestep: Time) -> Result<(),ThermalHydraulicsLibError> {

        self.pipe_fluid_array.advance_timestep_mut_self(timestep)?;
        self.pipe_shell.advance_timestep_mut_self(timestep)?;
        Ok(())
        
    }

    /// advances timestep for each HeatTransferEntity within the 
    /// NonInsulatedPipe
    ///
    /// gives each pipe the parallel tube treatment
    #[inline]
    pub fn advance_timestep(&mut self, 
    timestep: Time) -> Result<(),ThermalHydraulicsLibError> {

        // first, we need to perform timestep advancement 
        // like for: 
        // self.pipe_fluid_array.advance_timestep_mut_self(timestep)?;
        //
        // however, we must factor in the 1/number of tubes for each tube
        
        let one_over_number_of_tubes: f64 = 1.0/(self.number_of_tubes as f64);

        // we shall need to clone and convert the pipe_fluid_array into 
        // an actual fluid array 

        let mut fluid_array_clone: FluidArray = 
            self.pipe_fluid_array.clone().try_into()?;

        let mass_flowrate_over_all_tubes = 
            fluid_array_clone.get_mass_flowrate();

        // once the fluid array clone is done, we can advance timestep
        // I'll first copy code for advancing fluid array timesteps

        // there will always be at least 2 nodes

        let number_of_nodes = fluid_array_clone.len();
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
        fluid_array_clone.back_single_cv.rate_enthalpy_change_vector.clone();

        let front_cv_rate_enthalpy_change_vector: Vec<Power> = 
        fluid_array_clone.front_single_cv.rate_enthalpy_change_vector.clone();


        // compute power source for back node,
        // then multiply by one_over_number_of_tubes

        let mut total_enthalpy_rate_change_back_node = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            back_cv_rate_enthalpy_change_vector.iter() {

                total_enthalpy_rate_change_back_node += *enthalpy_chg_rate;
            }

        total_enthalpy_rate_change_back_node *= one_over_number_of_tubes;

        // then the front node, do the same thing

        let mut total_enthalpy_rate_change_front_node = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            front_cv_rate_enthalpy_change_vector.iter() {

                total_enthalpy_rate_change_front_node += *enthalpy_chg_rate;
            }

        total_enthalpy_rate_change_front_node *= one_over_number_of_tubes;

        // this front and back nodes will be an extra term added to the 
        // heat source vector S
        //
        // We need a blank temperature array to start the iteration 
        // process, so we just copy the old temperature array over 
        // as the initial guess (i think?)
        //

        let new_temperature_array: Array1<ThermodynamicTemperature>;
        // now let's start calculation 
        //
        let mut coefficient_matrix: Array2<ThermalConductance> = 
        Array::zeros((number_of_nodes, number_of_nodes));

        let mut power_source_vector: 
        Array1<Power> = Array::zeros(number_of_nodes);

        // ascertain if we have forward flow (doesn't change if we divide 
        // by the number of tubes

        let forward_flow: bool = mass_flowrate_over_all_tubes.ge(
            &MassRate::zero());

        // obtain some important parameters for calculation
        let material = fluid_array_clone.material_control_volume;
        let pressure = fluid_array_clone.pressure_control_volume;
        let bulk_temperature = fluid_array_clone.try_get_bulk_temperature()?;
        let total_volume = fluid_array_clone.get_component_length() *  
            fluid_array_clone.get_cross_sectional_area_immutable();
        let dt = timestep;
        let node_length = fluid_array_clone.get_component_length() 
            / number_of_nodes as f64;


        // for the fluid array, we do not keep track of node enthalpies,
        // instead, we keep track of temperatures and then calculate 
        // enthalpy changes using rho_cp calculated at the current 
        // timestep
        //
        // the m_cp must in turn be divided by the number of parallel tubes
        let volume_fraction_array: Array1<f64> = 
        fluid_array_clone.volume_fraction_array.iter().map(
            |&vol_frac| {
                vol_frac
            }
        ).collect();

        // rho_cp is determined by the temperature array 
        // m_cp will be determined by number of tubes (later)
        let rho_cp: Array1<VolumetricHeatCapacity> = 
        fluid_array_clone.temperature_array_current_timestep.iter().map(
            |&temperature| {
                try_get_rho_cp(material, temperature, pressure).unwrap()
            }
        ).collect();
        // energy balance for single tube is: 
        // m c_p dT/dt = -\sum H (T - T_lateral) - m_flow h_fluid(T) 
        // + m_flow h_fluid(adjacent T) + q
        //
        // of all these terms, only the m cp dT/dt term and HT are are 
        // implicitly calculated
        // 
        // We separate this out to get:
        //
        // m cp T / dt + \sum HT = 
        //
        // \sum HT_lateral - m_flow h_fluid(T_old) 
        // + m_flow h_fluid(adjacent T_old) + m cp / dt (Told)
        // + q
        //
        // so we will need to determine sum H and sum HT_lateral
        // as sum H is the relevant coefficient in the coefficient_matrix
        //
        //
        // When parallel


        //self.pipe_shell.advance_timestep_mut_self(timestep)?;
        Ok(())
        
    }


    /// advances timestep by spawning a thread 
    /// 
    pub fn advance_timestep_thread_spawn(&self,
        timestep: Time,) -> JoinHandle<Self> {

        // make a clone
        let mut heater_clone = self.clone();

        // move ptr into a new thread 

        let join_handle = thread::spawn(
            move || -> Self {


                // carry out the connection calculations
                heater_clone.advance_timestep(timestep).unwrap();
                
                heater_clone

            }
        );

        return join_handle;

    }
}
