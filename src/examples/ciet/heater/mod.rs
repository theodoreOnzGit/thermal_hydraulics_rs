//! The heater module here represents components from 
//!
//! BT-11 to BT-12
//! This is because heater inlet and outlet temperatures are measured 
//! using BT-11 and BT-12
//! 
//! However, BT-11 is located before the heater bottom head 
//! and BT-12 is located after the static mixer MX-10 
//!
//! Hence, we have four sections,
//!
//! the heater top head, heater bottom head, heated section,
//! mixer MX-10 and the static mixer pipe attached to MX-10 modelled 
//! in the Relap and SAM model
//!
//! So there is not only some residence time, but also other mechanisms 
//! for parasitic heat loss 
//!
//! Dr Dane De Wet's transform model does callibrate for this using 
//! a heat transfer coefficient of 20 W/(m^2 K) instead of 6 W/(m^2 K)
//!
//!
//! I intend to connect structural supports to the heater top and bottom 
//! head and callibrate the length of those structural supports 
//! as part of model callibration such that the heater inlet is 
//! 80C approximately, and the heater outlet is 102.45 C 
//!
//! at nominal heater power of 8 kW
//!
//! For this, I also want to ensure that the code runs fast enough,
//! at least faster than real time, so it is suitable for digital 
//! twin applications
//!
//!
//! 
//!
/// represents heater version 1, 
///
/// it allows fluid to flow through it in an annular tube
pub struct HeaterVersion1;

/// represents heater version 2, 
///
/// it has twisted tape
pub struct HeaterVersion2;

pub mod heater_version_2_bare;
pub use heater_version_2_bare::*;

pub mod heater_top_and_bottom_head_bare;
pub use heater_top_and_bottom_head_bare::*;

pub mod static_mixer_mx_10;
pub use static_mixer_mx_10::*;

use thermal_hydraulics_rs::prelude::alpha_nightly::*;
use uom::{si::{time::second, power::kilowatt}, ConstZero};

pub fn example_heater(){

    // construct structs

    let _heater_v1 = HeaterVersion1{};
    let _heater_v2 = HeaterVersion2{};


    // bare heater example
    let initial_temperature: ThermodynamicTemperature = 
    ThermodynamicTemperature::new::<degree_celsius>(79.12);
    let inlet_temperature = initial_temperature;
    let ambient_air_temp: ThermodynamicTemperature = 
    ThermodynamicTemperature::new::<degree_celsius>(21.76);

    let number_of_temperature_nodes: usize = 8;
    
    let mut heater_v2_bare = HeaterVersion2Bare::new_dewet_model(
        initial_temperature,
        ambient_air_temp,
        number_of_temperature_nodes
    );

    let mut inlet_bc: HeatTransferEntity = BCType::new_const_temperature( 
        inlet_temperature).into();

    let mut outlet_bc: HeatTransferEntity = BCType::new_adiabatic_bc().into();

    // time settings 

    let max_time = Time::new::<second>(75.0);
    let timestep = Time::new::<second>(0.01);
    let mut simulation_time = Time::ZERO;
    let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
    let heater_power = Power::new::<kilowatt>(8.0);

    // main loop
    
    while max_time > simulation_time {

        // create interactions 

        // inlet fluid density 

        let inlet_fluid_density: MassDensity = 
        LiquidMaterial::TherminolVP1.density(
            inlet_temperature).unwrap();

        // first node of heater fluid density 

        let therminol_array_clone: FluidArray 
        = heater_v2_bare.therminol_array.clone().try_into().unwrap();

        let therminol_array_temperature: Vec<ThermodynamicTemperature> = 
        therminol_array_clone.get_temperature_vector().unwrap();

        let steel_array_clone: SolidColumn 
        = heater_v2_bare.steel_shell.clone().try_into().unwrap();

        let steel_array_temperature: Vec<ThermodynamicTemperature> = 
        steel_array_clone.get_temperature_vector().unwrap();

        let back_cv_temperature: ThermodynamicTemperature = 
            therminol_array_temperature[0];

        let front_cv_temperature: ThermodynamicTemperature = 
        *therminol_array_temperature.iter().last().unwrap();

        let back_cv_density: MassDensity = 
        LiquidMaterial::TherminolVP1.density(
            back_cv_temperature).unwrap();

        let front_cv_density: MassDensity = 
        LiquidMaterial::TherminolVP1.density(
            front_cv_temperature).unwrap();

        // probably want to make this bit a little more user friendly
        let inlet_interaction: HeatTransferInteractionType = 
        HeatTransferInteractionType::new_advection_interaction(
            mass_flowrate,
            inlet_fluid_density,
            back_cv_density);

        let outlet_interaction = 
        HeatTransferInteractionType::new_advection_interaction(
            mass_flowrate,
            front_cv_density,
            front_cv_density,
        );

        // make axial connections to BCs 

        heater_v2_bare.therminol_array.link_to_back(
            &mut inlet_bc,
            inlet_interaction
        ).unwrap();

        heater_v2_bare.therminol_array.link_to_front(
            &mut outlet_bc,
            outlet_interaction
        ).unwrap();

        // make other connections 

        heater_v2_bare.lateral_and_miscellaneous_connections(
            mass_flowrate,
            heater_power
        );

        // calculate timestep
        heater_v2_bare.advance_timestep(
            timestep);


        simulation_time += timestep;

        // print outlet temperature 
        println!("Exit Temp {:?}", front_cv_temperature);

        // print surface temperature 
        println!("Steel array Temp: \n {:?} \n", steel_array_temperature);

        // print surface temperature 
        println!("Therminol Array Temp: \n{:?}", therminol_array_temperature);
    }

    // once simulation completed, write data


    todo!()

}
#[test]
pub fn example_heater_test(){

}
