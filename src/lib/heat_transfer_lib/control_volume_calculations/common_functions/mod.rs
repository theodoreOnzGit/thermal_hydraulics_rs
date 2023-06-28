/// thermal resistance models which are used for 
/// 1D steady state heat conduction problems.
pub mod thermal_resistance;
pub use thermal_resistance::*;


/// courant number, we want to help calculate stability over here
pub mod courant_number;
pub use courant_number::*;


extern crate uom;
use uom::si::f64::*;
/// This function calculates the formula:
/// H_cv (t+1)  = dt * (H_in - H_out + Q_s + W_s) + H_cv (t)
///
/// It is up to user discretion whether to use the enthalpy
/// in and out of the system and work done on the system
/// as well as heat supplied at current timestep or next timestep
///
#[inline]
pub fn get_control_volume_enthalpy_next_timestep(
    timestep: Time,
    enthalpy_out: Power,
    enthalpy_in: Power,
    heat_supplied_to_system: Power,
    work_done_on_system: Power,
    control_volume_enthalpy_current_timestep: Energy,
    ) -> Energy {

    let control_volume_enthalpy_next_timestep: Energy 
        = timestep * (enthalpy_in - 
                      enthalpy_out +
                      heat_supplied_to_system +
                      work_done_on_system) 
        + control_volume_enthalpy_current_timestep;

    return control_volume_enthalpy_next_timestep;

}

/// calculates enthalpy flow from
/// mass flowrate and specific energy (or available energy)
///
/// I will need to also calculate specific energy given a temperature.
/// This is fluid specific and we will need to code
/// correlations to convert specific energy (or available energy)
/// to temperature and vice versa
pub fn calculate_enthalpy_flow(
    m : MassRate,
    h : AvailableEnergy) -> Power {

    // m is mass flowrate
    // h is enthalpy per unit mass
    //
    // we are returning m * h which is enthalpy flowrate

    return m * h ;
}


/// a trait for temperature and specific
/// enthalpy conversion, not structly necessary to use
pub trait TemperatureEnthalpyConversion {

    /// obtains specific enthalpy given a fluid temperature
    fn temperature_to_enthalpy(
        fluid_temp: ThermodynamicTemperature) -> AvailableEnergy;

    /// obtains fluid temperature given a specific enthalpy
    fn enthalpy_to_temperature(
        fluid_enthalpy: AvailableEnergy) -> ThermodynamicTemperature;
}

/// calculate convection heat flux power input
///
/// Q = h (T_surface - T_fluid) A
pub fn calculate_convection_heat_flux_power_input(
    heat_transfer_coefficient : HeatTransfer,
    surface_temperature : ThermodynamicTemperature,
    fluid_temperature : ThermodynamicTemperature,
    surface_area : Area) -> Power {

    // note, i do this to calculate
    // delta T = surface_temperature - fluid_temperature
    let temperature_diff = 
        subtract_two_thermodynamic_temperatures(
        surface_temperature, fluid_temperature);

    return heat_transfer_coefficient * (temperature_diff) * surface_area; 

}





/// calculate overall heat flux power input
///
/// assuming a fixed surrounding temperature
///
/// calculates heat INPUT into fluid based on surrounding temperature
/// estimated fluid inlet and fluid outlet temperature
///
/// Q = U * A * LMTD
///
/// LMTD = (delta T in - delta T out) / (ln delta T in - ln delta T out)
///
/// U is overall_heat_transfer_coeff
///
/// A is the surface_area 
/// The surface_area you use can be the surface area of the inner 
/// or outer region of the pipe. BUT, the overall_heat_transfer_coeff 
/// must be adjusted accordingly
///
///
pub fn calculate_overall_heat_flux_power_input(
    overall_heat_transfer_coeff : HeatTransfer,
    ambient_temperature : ThermodynamicTemperature,
    fluid_temperature_in : ThermodynamicTemperature,
    fluid_temperature_out: ThermodynamicTemperature,
    surface_area : Area) -> Power {

    if fluid_temperature_in.value == fluid_temperature_out.value {
        return overall_heat_transfer_coeff * surface_area * 
            subtract_two_thermodynamic_temperatures(
            ambient_temperature , fluid_temperature_in);
    }
    // note, i do this to calculate
    // delta T = ambient_temperature - fluid_temperature
    let log_mean_temp_diff = 
        log_mean_temperature_difference(
            ambient_temperature,
            ambient_temperature, 
            fluid_temperature_in,
            fluid_temperature_out);

    return overall_heat_transfer_coeff * (log_mean_temp_diff) * surface_area; 

}



use uom::si::{temperature_interval, thermodynamic_temperature};


/// this function is a handy tool to convert temperature
/// interval types to thermodynamic temperatures
///
/// So whatever the temperature interval given, I will 
/// add 273.15K to it to make it a proper thermodynamic
/// Temperature
///
/// For reference I will be testing that 83F is 28.333C
///
/// However, an 83F temperature interval is actually 
/// 46.11 K.
///
/// I need to subtract 32F off first so that the temperature
/// interval of 51F represents 28.33 C
///
/// ```rust
/// extern crate approx;
/// 
/// use thermal_hydraulics_rs::heat_transfer_lib
/// ::control_volume_calculations::common_functions;
///
/// use uom::si::{temperature_interval, thermodynamic_temperature};
/// use uom::si::f64::*;
///
/// let ti = TemperatureInterval::new::
/// <temperature_interval::degree_fahrenheit>(83_f64-32_f64);
///
///
/// let expected_temp = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_celsius>(28.333333);
///
/// let test_temp = common_functions::
/// temperature_interval_to_thermodynamic_temperature
/// (ti);
///
/// approx::assert_relative_eq!(expected_temp.value, test_temp.value, 
/// max_relative=0.001);
///
/// ```
pub fn temperature_interval_to_thermodynamic_temperature(
    temp: TemperatureInterval) -> ThermodynamicTemperature {

    // first we obtain the temperature interval in kelvin
    // usually values in the Quantity Structs are stored as kelvin
    // we'll need to demo this
    let temperature_interval_value : f64 = 
        temp.value;

    // second we add 273.15K to the interval to convert it
    // to a temperature

    let temperature_absolute_value : f64 = temperature_interval_value +
        273.15;

    // third we return this value as a dimensioned thermodynamic
    // temperature

    let temp_return = ThermodynamicTemperature::new::
        <thermodynamic_temperature::kelvin>(temperature_absolute_value);

    return temp_return;

}


/// Converts a thermodynamic temperature value to a 
/// temperature interval by subtracting 273.15K from the value
///
/// For reference I will be testing that 83F is 28.333C
///
/// However, if i want 83F eventually, i need to supply a 
/// temperature interval of 83 F - 32 F =  51 F
///
/// 51 F + 273.15K = 28.333C
/// Of course, you'll need to do all the unit conversions and etc.
///
/// ```rust
/// extern crate approx;
/// 
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// control_volume_calculations::common_functions;
///
/// use uom::si::{temperature_interval, thermodynamic_temperature};
/// use uom::si::f64::*;
///
/// let tt = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_fahrenheit>(83_f64);
///
///
/// let expected_temp = TemperatureInterval::new::
/// <temperature_interval::degree_celsius>(28.333333);
///
/// let test_temp = common_functions::
/// thermodynamic_temperature_to_temperature_interval
/// (tt);
///
/// approx::assert_relative_eq!(expected_temp.value, test_temp.value, 
/// max_relative=0.001);
///
/// ```
pub fn thermodynamic_temperature_to_temperature_interval(
    temp: ThermodynamicTemperature) -> TemperatureInterval {

    // first we obtain the thermodynamic temperature in kelvin
    // usually values in the Quantity Structs are stored as kelvin
    // we'll need to demo this
    let temperature_absolute_value : f64 = 
        temp.value;

    // second we subtract 273.15K from the absolute temp to convert it
    // to a temperature interval

    let temperature_interval_value : f64 = temperature_absolute_value -
        273.15;

    // third we return this value as a dimensioned 
    // temperature interval

    let temp_return = TemperatureInterval::new::
        <temperature_interval::kelvin>(temperature_interval_value);

    return temp_return;
}


/// subtracts two thermodynamic temperatures from each
/// other to obtain a temperature interval
///
/// two values are supplied, t1 and t2
///
/// So i'm going to subtract 83F from 600K
/// 83 F is 301.5 K approximately
///
///
/// ```rust
/// extern crate approx;
/// use thermal_hydraulics_rs::heat_transfer_lib
/// ::control_volume_calculations::common_functions;
///
/// use uom::si::{temperature_interval, thermodynamic_temperature};
/// use uom::si::f64::*;
///
/// let t1 = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::kelvin>(600_f64);
///
/// let t2 = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_fahrenheit>(83_f64);
///
/// let expected_temp_value = t1.value - t2.value;
///
/// let test_temp = common_functions::subtract_two_thermodynamic_temperatures(
/// t1,t2);
/// approx::assert_relative_eq!(expected_temp_value, test_temp.value, 
/// max_relative=0.001);
///
/// ```
///
///
///
pub fn subtract_two_thermodynamic_temperatures(
    t1: ThermodynamicTemperature,
    t2: ThermodynamicTemperature) -> TemperatureInterval {

    let temperature_interval_value = t1.value - t2.value;

    return TemperatureInterval::new::<temperature_interval::
        kelvin>(temperature_interval_value);

}


/// LMTD = (delta T in - delta T out) / (ln delta T in - ln delta T out)
///
/// note that reversing the order of delta T in and out doesn't really
/// matter, as long as both numerator and denominator are reversed
/// correctly
///
/// However, hot fluid temperatures and cold fluid temperature CANNOT
/// be mixed up, otherwise the logarithms will return an error
///
/// ```rust
/// extern crate approx;
/// use thermal_hydraulics_rs::heat_transfer_lib
/// ::control_volume_calculations::common_functions;
///
/// use uom::si::{temperature_interval, thermodynamic_temperature};
/// use uom::si::f64::*;
///
/// let cold_fluid_temp_A = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_celsius>(21.0);
///
/// let cold_fluid_temp_B = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_celsius>(20.0);
///
/// let hot_fluid_temp_A = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_celsius>(48.0);
///
/// let hot_fluid_temp_B = ThermodynamicTemperature::new::
/// <thermodynamic_temperature::degree_celsius>(50.0);
///
/// let A_temperature_interval_value : f64 = hot_fluid_temp_A.value - 
/// cold_fluid_temp_A.value;
///
/// let B_temperature_interval_value : f64 = hot_fluid_temp_B.value - 
/// cold_fluid_temp_B.value;
///
/// let mut LMTD_value_expected = 
/// (A_temperature_interval_value - B_temperature_interval_value)/
/// (A_temperature_interval_value.ln() - 
/// B_temperature_interval_value.ln());
///
/// let LMTD_test = common_functions::log_mean_temperature_difference(
/// cold_fluid_temp_A,
/// cold_fluid_temp_B,
/// hot_fluid_temp_A,
/// hot_fluid_temp_B);
///
///
/// approx::assert_relative_eq!(LMTD_value_expected, LMTD_test.value, 
/// max_relative=0.001);
///
/// // test 2 makes it more obvious
///
/// let mut LMTD_value_expected = 
/// ((48_f64 - 21_f64) - (50_f64-20.0))/
/// ((48_f64-21_f64).ln() - 
/// (50_f64-20.0).ln());
///
/// approx::assert_relative_eq!(LMTD_value_expected, LMTD_test.value, 
/// max_relative=0.001);
///
/// ```
pub fn log_mean_temperature_difference(
    temp_cold_fluid_a: ThermodynamicTemperature,
    temp_cold_fluid_b: ThermodynamicTemperature,
    temp_hot_fluid_a: ThermodynamicTemperature,
    temp_hot_fluid_b: ThermodynamicTemperature) -> TemperatureInterval {

    if temp_hot_fluid_a.value < temp_cold_fluid_a.value {
        panic!("LMTD: hot fluid temperature input colder than \n
               cold fluid temperature input")
    }

    if temp_hot_fluid_b.value < temp_cold_fluid_b.value {
        panic!("LMTD: hot fluid temperature input colder than \n
               cold fluid temperature input")
    }
    
    let a_temperature_interval = 
        subtract_two_thermodynamic_temperatures(
            temp_hot_fluid_a, temp_cold_fluid_a);

    let b_temperature_interval = 
        subtract_two_thermodynamic_temperatures(
            temp_hot_fluid_b, temp_cold_fluid_b);

    let numerator = a_temperature_interval - b_temperature_interval;

    let denominator = a_temperature_interval.value.ln() - 
        b_temperature_interval.value.ln();


    return numerator/denominator;

}

