use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;
use ndarray_linalg::{Solve, error};
use uom::si::f64::*;
use uom::ConstZero;
use ndarray::*;
use uom::si::thermodynamic_temperature::kelvin;


impl ArrayCVType {

    /// sets the current temperature vector and other 
    /// properties based on temperatures calculated
    /// for the next timestep
    /// and also other cleanup work
    pub fn advance_timestep(&mut self, timestep: Time) -> Result<(),String>{

        match self {
            ArrayCVType::Cartesian1D(cartesian_1d_cv) => {
                cartesian_1d_cv.advance_timestep(timestep)
            },
        }
    }
}
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn solve_conductance_matrix_power_vector(
    thermal_conductance_matrix: Array2<ThermalConductance>,
    power_vector: Array1<Power>)
-> Result<Array1<ThermodynamicTemperature>, error::LinalgError>{

    // I can of course convert it into f64 types 
    //
    //

    let get_value_conductance = |conductance: &ThermalConductance| {
        return conductance.value;
    };

    let get_value_power = |power: &Power| {
        return power.value;
    };

    // i'm allowing non snake case so that the syntax is the same as 
    // GeN-Foam
    #[allow(non_snake_case)]
    let M: Array2<f64> = 
    thermal_conductance_matrix.map(get_value_conductance);

    #[allow(non_snake_case)]
    let S: Array1<f64> = power_vector.map(get_value_power);

    // now for the raw temperature matrix 

    #[allow(non_snake_case)]
    let T: Array1<f64> = M.solve(&S)?;

    // To check for unit safety, I can just perform one calc

    let _unit_check: Power = 
    power_vector[0] + 
    thermal_conductance_matrix[[0,0]] 
    * ThermodynamicTemperature::ZERO;

    // now map T back to a ThermodynamicTemperature
    // T is already a ThermodynamicTemperature, so don't need to manually 
    // convert, do it in kelvin 
    //

    let convert_f64_to_kelvin_temperature = |float: &f64| {
        return ThermodynamicTemperature::new::<kelvin>(*float);
    };


    // this is the last step
    let temperature_vector: Array1<ThermodynamicTemperature> 
    = T.map(convert_f64_to_kelvin_temperature);

    return Ok(temperature_vector);
}
