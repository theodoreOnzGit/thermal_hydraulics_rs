use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::specific_enthalpy::try_get_temperature_from_h;

use super::SingleCVNode;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;

impl SingleCVNode {
    /// this function performs necessary calculations to move 
    /// the state of the control volume to the next time step
    ///
    /// calculates the new enthalpy of the 
    /// and cleans out the all power vectors and time step vectors
    #[inline]
    pub fn advance_timestep(&mut self, timestep: Time) -> Result<(), ThermalHydraulicsLibError>{


        // first thing is to sum up all enthalpy changes
        let mut total_enthalpy_rate_change = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            self.rate_enthalpy_change_vector.iter() {

                total_enthalpy_rate_change += *enthalpy_chg_rate;
            }

        let enthalpy_next_timestep = total_enthalpy_rate_change * 
        timestep +
        self.current_timestep_control_volume_specific_enthalpy
            * self.mass_control_volume;

        let specific_enthalpy_next_timestep = 
        enthalpy_next_timestep/self.mass_control_volume;


        self.next_timestep_specific_enthalpy 
            = specific_enthalpy_next_timestep;

        // at the end of each timestep, set 
        // current_timestep_control_volume_specific_enthalpy
        // to that of the next timestep

        self.current_timestep_control_volume_specific_enthalpy
            = specific_enthalpy_next_timestep;

        // now things get a bit clunky here, I will need to set 
        // the mass of the control volume by the temperature after 
        // I calculated it using control volume mass from the last 
        // timestep. For large changes in density, this may not 
        // work well (instabilities??) but for liquids, I hope it's 
        // okay

        self.set_liquid_cv_mass_from_temperature()?;
        // clear the enthalpy change vector and timestep vector 
        // also the mass flowrate vector

        self.clear_vectors().unwrap();
        // increase timestep (last step)

        // set temperatures for the cv

        let _ = self.get_temperature_from_enthalpy_and_set()?;

        return Ok(());
    }
    /// clears all vectors for next timestep
    /// This is important for the advance timestep method
    pub fn clear_vectors(&mut self) 
    -> Result<(), ThermalHydraulicsLibError>{


        self.rate_enthalpy_change_vector.clear();
        self.max_timestep_vector.clear();
        self.volumetric_flowrate_vector.clear();

        Ok(())
    }


    /// this function calculates the conductance interaction between this 
    /// control volume and another control volume
    ///
    /// In the case of advection, the other control volume is placed in front 
    /// of this control volume
    #[inline]
    pub fn calculate_conductance_interaction_to_front_singular_cv_node(
        &mut self,
        single_cv_2: &mut SingleCVNode,
        interaction: HeatTransferInteractionType)-> Result<(), ThermalHydraulicsLibError>{

        // let's get the two temperatures of the control volumes first
        // so let me get the enthalpies, and then their respective 
        // temperatures 

        let single_cv_1_enthalpy = self.
            current_timestep_control_volume_specific_enthalpy;
        let single_cv_2_enthalpy = single_cv_2.
            current_timestep_control_volume_specific_enthalpy;

        // to get the temperatures, we'll need the material as well 
        let single_cv_1_material = self.material_control_volume;
        let single_cv_2_material = single_cv_2.material_control_volume;

        // we'll also need to get their pressures 
        let single_cv_1_pressure = self.pressure_control_volume;
        let single_cv_2_pressure = single_cv_2.pressure_control_volume;

        // we will now get their respective temperatures 
        //
        // (note, this is extremely computationally expensive as it 
        // is iterative in nature)
        //
        // two solutions here, 
        // one: store cv temperature in single cv, 
        // so it can be readily accessed
        //
        // two: cheaper method of getting t from h.


        let single_cv_1_temperature: ThermodynamicTemperature;
        let single_cv_2_temperature: ThermodynamicTemperature;

        let experimental_code = true; 
        if experimental_code {

            single_cv_1_temperature = self.temperature;
            single_cv_2_temperature = single_cv_2.temperature;

        } else {

            // original code
            single_cv_1_temperature = try_get_temperature_from_h(
                single_cv_1_material, 
                single_cv_1_enthalpy, 
                single_cv_1_pressure)?;
            single_cv_2_temperature = try_get_temperature_from_h(
                single_cv_2_material, 
                single_cv_2_enthalpy, 
                single_cv_2_pressure)?;

        }

        // now that we got their respective temperatures we can calculate 
        // the thermal conductance between them
        //
        // for conduction for instance, q = kA dT/dx 
        // conductance is watts per kelvin or 
        // q = (kA)/dx * dT
        // conductance here is kA/dx
        // thermal resistance is 1/conductance
        //
        // for convection, we get: 
        // q = hA (Delta T)
        // hA becomes the thermal conductance
        //
        // If we denote thermal conductance as Htc
        // 
        // Then a general formula for heat flowing from 
        // temperature T_1 to T_2 is 
        //
        // T_1 --> q --> T_2 
        //
        // q = - Htc (T_2 - T_1)

        // 
        let thermal_conductance = interaction.get_thermal_conductance_based_on_interaction(
            single_cv_1_temperature, 
            single_cv_2_temperature,
            single_cv_1_pressure, 
            single_cv_2_pressure)?;

        // suppose now we have thermal conductance, we can now obtain the 
        // power flow
        //

        let cv_2_temp_minus_cv_1_temp_kelvin: f64 = 
            single_cv_2_temperature.get::<kelvin>() - 
            single_cv_1_temperature.get::<kelvin>();

        let cv_2_temp_minus_cv_1: TemperatureInterval = 
            TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(
                cv_2_temp_minus_cv_1_temp_kelvin);

        let heat_flowrate_from_cv_1_to_cv_2: Power = 
            - thermal_conductance * cv_2_temp_minus_cv_1;

        // now, we add a heat loss term to cv_1 
        // and a heat gain term to cv_2 
        //
        // using timestep
        // the signs should cancel out

        self.rate_enthalpy_change_vector.
            push(-heat_flowrate_from_cv_1_to_cv_2);
        single_cv_2.rate_enthalpy_change_vector.
            push(heat_flowrate_from_cv_1_to_cv_2);


        // for solids mesh fourier number need only 
        // be done once, not every time 
        // an interaction is formed 
        //
        // probably the cell stability fourier number will be done in the 
        // constructor. however, with convection, the time scale must be 
        // recalculated at every time step. so it really depends whether 
        // it's solid or fluid control volume
        // Actually, for solid CV, I will also need to recalculate time scale 
        // based on the material thermal thermal_diffusivity
        //
        // For liquid CV, still need to calculate time scale based 
        // on convection flow
        // match statement is meant to tell that liquid CVs are not quite 
        // ready for use
        //
        // but the liquid timescales are calculated at the cv level only 
        // after all volumetric flowrates are calculated
        //
        // so don't really need any new timescale calculations

        return Ok(());

    }
}

