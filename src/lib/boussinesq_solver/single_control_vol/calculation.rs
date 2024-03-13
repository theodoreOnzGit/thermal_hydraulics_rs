use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::advection_heat_rate;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::boussinesq_solver::heat_transfer_correlations::
heat_transfer_interactions::heat_transfer_interaction_enums::DataAdvection;
use crate::boussinesq_solver::heat_transfer_correlations::
heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::
specific_enthalpy::try_get_temperature_from_h;

use super::SingleCVNode;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;
use uom::num_traits::Zero;

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

}

