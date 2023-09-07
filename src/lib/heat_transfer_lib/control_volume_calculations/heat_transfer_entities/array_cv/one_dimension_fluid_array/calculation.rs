use ndarray::*;
use ndarray_linalg::error::LinalgError;
use uom::si::f64::*;
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

    /// clears all vectors for next timestep
    /// This is important for the advance timestep method
    pub fn clear_vectors(&mut self) 
        -> Result<(), ThermalHydraulicsLibError>{

        self.lateral_adjacent_array_conductance_vector.clear();
        self.lateral_adjacent_array_temperature_vector.clear();
        Ok(())
    }
}
