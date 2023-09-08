use super::SolidColumn;
use uom::si::f64::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
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

        // now that we've gotten all the important properties, we can 
        // start matrix construction

        
        
        todo!()
    }
}
