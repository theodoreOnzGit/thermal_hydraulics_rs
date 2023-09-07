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

    /// connects an adjacent solid or fluid node laterally 
    /// with a given average thermal conductance
    /// note that doing so with 
    pub fn lateral_link_new_temperature_vector_avg_conductance(&mut self,
    average_thermal_conductance: ThermalConductance,
    temperature_vec: Vec<ThermodynamicTemperature>) 
        -> Result<(), ThermalHydraulicsLibError>{

        let number_of_temperature_nodes = NUMBER_OF_NODES;

        // check if temperature_vec has the correct number_of_temperature_nodes

        if temperature_vec.len() !=  number_of_temperature_nodes {
            let shape_error = ShapeError::from_kind(
                ErrorKind::IncompatibleShape
            );

            let linalg_error = LinalgError::Shape(shape_error);

            return Err(ThermalHydraulicsLibError::LinalgError
                (linalg_error));

        }

        // now let's make a new temperature array 

        let mut temperature_arr: [ThermodynamicTemperature; NUMBER_OF_NODES]
        = [ThermodynamicTemperature::new::<kelvin>(0.0); NUMBER_OF_NODES];

        // assign the temperatures 
        
        for (idx, temperature) in temperature_vec.iter().enumerate() {

            temperature_arr[idx] = *temperature;

        }

        // push it to the lateral adjacent array temp vec 

        self.lateral_adjacent_array_temperature_vector.push(temperature_arr);

        // next, construct conductance vector

        let conductance_arr: [ThermalConductance; NUMBER_OF_NODES]
        = [average_thermal_conductance; NUMBER_OF_NODES];

        self.lateral_adjacent_array_conductance_vector.push(conductance_arr);

        Ok(())
    }

}
