use std::f64::consts::PI;
use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;

use peroxide::prelude::erfc;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::data_enum_structs::DataUserSpecifiedConvectionResistance;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::{link_heat_transfer_entity, HeatTransferInteractionType};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::SolidMaterial::{SteelSS304L, Copper};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::Material;
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, OuterDiameterThermalConduction, SurfaceArea, SingleCVNode, CVType};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::CVType::SingleCV;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::density::density;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy};



use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::thermal_diffusivity::thermal_diffusivity;
use uom::si::f64::*;
use uom::si::length::{centimeter, meter};
use uom::si::power::watt;
use uom::si::ratio::ratio;
use uom::si::temperature_interval::degree_celsius as interval_deg_c;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

/// This is an example of transient conduction case where analytical 
/// solutions have been well know 
///
/// This is conduction in a semi infinite medium with constant 
/// temperature boundary conditions 
///
/// These analytical solutions are well known
/// and contain error functions 
///
/// For constant temperature, the temperature at point x 
/// and time t is given by 
///
/// theta (x,t) = erfc (x / (2.0 * sqrt{alpha t}) )
///
///
/// Trojan, M. "Transient Heat Conduction in Semiinfinite 
/// Solid with Surface Convection." Encyclopedia of (2014).
///
/// One may note that the (x / (2.0 * sqrt{alpha t}) )
/// is simply the reciprocal of a Fourier like number multiplied by 
/// a constant
///
/// Fo = (alpha t)/L^2 
///
/// In this case, the length scale is x and t is the time.
///
/// How long must the medium be for it to be semi infinite? 
///
/// erfc (zeta) can be plotted on a graph 
/// at erfc(zeta) where zeta >= 2.0,
/// erfc(zeta) <= 0.005 (0.5\% change)
///
/// theta (x,t) = erfc (1 / (2.0 * sqrt{Fo(x,t)}) )
///
///
/// probably use some published reference, better quality
/// // https://www.unipamplona.edu.co/unipamplona/portalIG/home_34/recursos/01general/21082014/unidad_2_termo_ii.pdf
///
///
#[test]
fn transient_conduction_semi_infinite_copper_medium() 
-> Result<(), String>{

    // let's first do the analytical solution


    let theta_error_fn = |fourier_number_x_t: Ratio| 
    -> Result<Ratio,String> {

        let fourier_value = fourier_number_x_t.value;

        // there's going to be a case where the Fourier number is exactly 
        // zero
        // We cannot divide things by zero 
        // but erfc of anything bigger than 2 is close to zero 
        // so if Fo -> 0, erfc(1/Fo) -> 0 
        // hence, a Fourier number of zero results in theta = 0 
        //

        if fourier_value == 0.0 {
            let theta_ratio = Ratio::new::<ratio>(0.0);
            return Ok(theta_ratio);
        }

        if fourier_value < 0.0 {

            return Err("negative fourier value".to_string());
        }


        // theta (x,t) = erfc (1 / (2.0 * sqrt{Fo(x,t)}) )
        let exponent_denominator = 2.0 * fourier_value.sqrt();

        let exponent: f64 = 1.0/exponent_denominator;

        let theta_value = erfc(exponent);

        let theta_ratio = Ratio::new::<ratio>(theta_value);

        return Ok(theta_ratio);

    };

    // let's do from t = 0 to t = 20 in 4 second intervals
    //
    // then we will print out the temperature profiles from x = 0 to 
    // x = 1m
    // Based on some preliminary calculations, the maximum length where 
    // 1/(2 sqrt Fo)  = 2 
    // at t = 20  
    // is about 0.2 m
    //
    // therefore, I'll place about 5 nodes there from x = 0m to 
    // x = 0.2m, these will record for us the temperature profile 
    // at a certain time

    let time_vector: Vec<Time> = vec![
        Time::new::<second>(0.0),
        Time::new::<second>(4.0),
        Time::new::<second>(8.0),
        Time::new::<second>(12.0),
        Time::new::<second>(16.0),
        Time::new::<second>(20.0),
    ];

    let length_vector: Vec<Length> = vec![
        Length::new::<meter>(0.0),
        Length::new::<meter>(0.05),
        Length::new::<meter>(0.10),
        Length::new::<meter>(0.15),
        Length::new::<meter>(0.20),
    ];

    // let's make the csv writer 

    use csv::Writer;
    let mut wtr = Writer::from_path("analytical_1d_transient_conduction.csv")
        .unwrap();

    wtr.write_record(&["time_seconds",
        "temp_0_0_meters", 
        "temp_0_5_meters",
        "temp_0_10_meters", 
        "temp_0_15_meters", 
        "temp_0_20_meters", 
    ]).unwrap();


    wtr.flush().unwrap();
    // also the copper thermal_diffusivity

    let copper = Material::Solid(Copper);
    let pressure = Pressure::new::<atmosphere>(1.0);
    let copper_initial_temperature = 
    ThermodynamicTemperature::new::<degree_celsius>(21.67);

    let boundary_condition_temperature = 
    ThermodynamicTemperature::new::<degree_celsius>(80.0);

    // note that diffusivity changes with temperature, but we shall not 
    // assume this is the case, and just obtain an approximate 
    // analytical solution 

    let copper_thermal_diffusivity_alpha: DiffusionCoefficient 
    = thermal_diffusivity(copper, copper_initial_temperature, pressure)?;

    for time in time_vector.iter() {

        // initialise a temperature vector 

        let mut temp_vector: Vec<ThermodynamicTemperature> = 
        vec![];

        // make a nested length loop 

        for length in length_vector.iter() {

            // if length is zero, then BC implies that temperature 
            // must be BC temperature 

            if length.value == 0.0 {
                temp_vector.push(boundary_condition_temperature);
            }

            // if time is zero, then initial conditions mean that 
            // temperature is the copper initial temperature

            else if time.value == 0.0 {
                temp_vector.push(copper_initial_temperature);
            } 
            else {

                // calc fourier number 
                let fourier_number: Ratio = 
                (copper_thermal_diffusivity_alpha * *time)
                / *length 
                / *length;

                let theta = theta_error_fn(fourier_number)?;

                // theta = (T(x,t) - T_i)/(T_BC - T_i)

                let temperature_diff: TemperatureInterval = 
                TemperatureInterval::new::<interval_deg_c>(
                    boundary_condition_temperature.value - 
                    copper_initial_temperature.value);

                let temperature_x_t: ThermodynamicTemperature = 
                copper_initial_temperature + 
                theta * temperature_diff;

                temp_vector.push(temperature_x_t);
            }

        }

        // once the temperature vector is finished, we can write it 
        // to csv file 

        let temperature_1 = temp_vector[0];
        let temperature_2 = temp_vector[1];
        let temperature_3 = temp_vector[2];
        let temperature_4 = temp_vector[3];
        let temperature_5 = temp_vector[4];

        wtr.write_record(&[&time.value.to_string(),
            & *temperature_1.value.to_string(),
            & *temperature_2.value.to_string(),
            & *temperature_3.value.to_string(),
            & *temperature_4.value.to_string(),
            & *temperature_5.value.to_string(),
        ]).unwrap();

        wtr.flush().unwrap();

        
    }

    // this setup is meant to be emulated using control volumes with 
    // some thermal resistances between them

    todo!("need to do 4 control vol and 1 BC for 1d transient conduction");
    Ok(())
}
