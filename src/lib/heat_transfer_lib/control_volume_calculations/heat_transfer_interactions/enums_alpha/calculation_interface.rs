use uom::si::f64::*;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::get_thermal_conductance_based_on_interaction;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::HeatTransferInteractionType;

impl HeatTransferInteractionType {

    /// attempts to get the thermal conductance from a heat transfer 
    /// interaction
    #[inline]
    pub fn try_get_thermal_conductance(&self,
    temperature_1: ThermodynamicTemperature,
    temperature_2: ThermodynamicTemperature,
    pressure_1: Pressure, 
    pressure_2: Pressure) -> Result<ThermalConductance
    ,ThermalHydraulicsLibError> {

        // todo: better error handling rather than unwrap
        let conductance: ThermalConductance = get_thermal_conductance_based_on_interaction(
            temperature_1,
            temperature_2,
            pressure_1,
            pressure_2,
            self.clone(),
        ).unwrap();

        Ok(conductance)
    }
}
