use super::SolidColumn;
use uom::si::f64::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use ndarray_linalg::error::LinalgError;
use ndarray::*;
use uom::si::power::watt;

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


        todo!()
    }
}
