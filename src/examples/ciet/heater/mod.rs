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
use core::time;
use std::thread::{self};
use std::thread::JoinHandle;
use std::time::SystemTime;

use csv::Writer;
pub use heater_version_2_bare::*;

pub mod heater_top_and_bottom_head_bare;
pub use heater_top_and_bottom_head_bare::*;

pub mod static_mixer_mx_10;
pub use static_mixer_mx_10::*;

pub mod struct_supports;

use thermal_hydraulics_rs::prelude::alpha_nightly::*;
use uom::ConstZero;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::length::centimeter;
use uom::si::length::foot;
use uom::si::power::kilowatt;
use uom::si::time::second;
use uom::si::time::minute;

use self::struct_supports::StructuralSupport;

pub mod heated_section_example;
pub mod heater_inclusive_top_bottom_head_example;
pub mod heater_plus_mx_10_without_supports;
pub mod heater_plus_mx_10_and_struct_support_example;

pub fn example_heater(){

    // construct structs

    let _heater_v1 = HeaterVersion1{};
    let _heater_v2 = HeaterVersion2{};


    // bare heater plus heads exaample
    let initial_temperature: ThermodynamicTemperature = 
    ThermodynamicTemperature::new::<degree_celsius>(79.12);
    let inlet_temperature = initial_temperature;
    let ambient_air_temp: ThermodynamicTemperature = 
    ThermodynamicTemperature::new::<degree_celsius>(21.67);

    let number_of_inner_temperature_nodes: usize = 6;
    
    let mut heater_v2_bare = HeaterVersion2Bare::new_dewet_model(
        initial_temperature,
        ambient_air_temp,
        number_of_inner_temperature_nodes
    );



    let mut heater_top_head_bare: HeaterTopBottomHead 
    = HeaterTopBottomHead::new_top_head(
        initial_temperature,
        ambient_air_temp);

    let mut heater_bottom_head_bare: HeaterTopBottomHead 
    = HeaterTopBottomHead::new_bottom_head(
        initial_temperature,
        ambient_air_temp);

    // calibration of heat transfer coeff
    let calibration_mode = true; 

    if calibration_mode {

        let h_to_air = HeatTransfer::new::<watt_per_square_meter_kelvin>
            (20.0);
        heater_v2_bare = HeaterVersion2Bare::_user_callibrated_htc_to_air_model(
            initial_temperature,
            ambient_air_temp,
            number_of_inner_temperature_nodes,
            h_to_air
        );

        heater_top_head_bare = HeaterTopBottomHead:: 
            _new_user_callibrated_top_head(
                initial_temperature,
                ambient_air_temp,
                h_to_air
            );
        heater_bottom_head_bare = HeaterTopBottomHead:: 
            _new_user_callibrated_bottom_head(
                initial_temperature,
                ambient_air_temp,
                h_to_air
            );
    }

    let mut static_mixer_mx_10_object: StaticMixerMX10 
    = StaticMixerMX10::new_static_mixer(
        initial_temperature,
        ambient_air_temp);

    let mut static_mixer_mx_10_pipe: StaticMixerMX10 
    = StaticMixerMX10::new_static_mixer_pipe(
        initial_temperature,
        ambient_air_temp);

    let struct_support_equiv_diameter: Length = Length::new::<inch>(0.5);
    let struc_support_equiv_length: Length = Length::new::<foot>(1.0);


    let mut structural_support_heater_top_head = 
    StructuralSupport::new_steel_support_cylinder(
        struc_support_equiv_length,
        struct_support_equiv_diameter,
        initial_temperature,
        ambient_air_temp);

    let mut structural_support_heater_bottom_head = 
    structural_support_heater_top_head.clone();

    let mut structural_support_mx_10 = 
    structural_support_heater_top_head.clone();

    let approx_support_conductance: ThermalConductance = 
    structural_support_heater_top_head.get_axial_node_to_bc_conductance();


    let support_conductance_interaction = HeatTransferInteractionType::
        UserSpecifiedThermalConductance(approx_support_conductance);


    let mut inlet_bc: HeatTransferEntity = BCType::new_const_temperature( 
        inlet_temperature).into();

    let mut outlet_bc: HeatTransferEntity = BCType::new_adiabatic_bc().into();

    let mut ambient_air_temp_bc: HeatTransferEntity = 
    inlet_bc.clone();

    // time settings 

    let max_time = Time::new::<minute>(2.0);
    // on my pc, the simulation time using 
    // cargo run --release 
    // is less than 10ms
    let timestep = Time::new::<second>(0.015);
    let mut simulation_time = Time::ZERO;
    let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
    let heater_power = Power::new::<kilowatt>(8.0);

    let loop_time = SystemTime::now();


    // csv writer
    let mut wtr = Writer::from_path("lib_heater_example_steady_state.csv")
        .unwrap();
    wtr.write_record(&["time_seconds",
        "heater_power_kilowatts",
        "therminol_temperature_celsius",
        "shell_temperature_celsius",
        "timestep_seconds",])
        .unwrap();

    let mut temp_profile_wtr = Writer::from_path(
        "lib_heater_example_temp_profile.csv")
        .unwrap();

    let number_of_nodes_heated_section = number_of_inner_temperature_nodes + 2;
    let node_length_heated_section = heater_v2_bare.get_component_length()
    / number_of_nodes_heated_section as f64;

    // for temperature profile
    let mut header_vec: Vec<String> = vec![];
    for index in 0..number_of_nodes_heated_section {

        let half_node_length = 0.5 * node_length_heated_section;
        let mid_node_length: Length = 
        index as f64 * node_length_heated_section + half_node_length;

        let prefix: String = "heater_temp_celsius_at_".to_string();

        let suffix: String = "_cm".to_string();

        let mid_node_length_cm: f64 = 
        mid_node_length.get::<centimeter>();

        let mid_node_length_string: String = 
        mid_node_length_cm.to_string();

        let header: String = prefix + &mid_node_length_string + &suffix;

        header_vec.push(header);


    }
    temp_profile_wtr.write_record(&header_vec).unwrap();

    // main loop
    
    let main_loop = thread::spawn( move || {
        while max_time > simulation_time {

            // time start 
            let loop_time_start = loop_time.elapsed().unwrap();
            // create interactions 


            // let's get heater temperatures for post processing
            // as well as the interaction
            // for simplicity, i use the boussineseq approximation,
            // which assumes that heat transfer is governed by 
            // average density (which doesn't change much for liquid 
            // anyway)


            let mut therminol_array_clone: FluidArray 
            = heater_v2_bare.therminol_array.clone().try_into().unwrap();

            let _therminol_array_temperature: Vec<ThermodynamicTemperature> = 
            therminol_array_clone.get_temperature_vector().unwrap();

            let heater_surface_array_clone: SolidColumn 
            = heater_v2_bare.steel_shell.clone().try_into().unwrap();

            let heater_surface_array_temp: Vec<ThermodynamicTemperature> = 
            heater_surface_array_clone.get_temperature_vector().unwrap();

            let heater_fluid_bulk_temp: ThermodynamicTemperature = 
            therminol_array_clone.try_get_bulk_temperature().unwrap();

            let heater_top_head_bare_therminol_clone: FluidArray = 
            heater_top_head_bare.therminol_array.clone().try_into().unwrap();

            let heater_top_head_exit_temperature: ThermodynamicTemperature = 
            heater_top_head_bare_therminol_clone.get_temperature_vector()
                .unwrap().into_iter().last().unwrap();

            let static_mixer_therminol_clone: FluidArray = 
            static_mixer_mx_10_object.therminol_array.clone().try_into().unwrap();

            let _static_mixer_exit_temperature: ThermodynamicTemperature
            = static_mixer_therminol_clone.get_temperature_vector().unwrap()
                .into_iter().last().unwrap();

            let static_mixer_pipe_therminol_clone: FluidArray = 
            static_mixer_mx_10_pipe.therminol_array.clone().try_into().unwrap();


            let heater_therminol_avg_density: MassDensity = 
            LiquidMaterial::TherminolVP1.density(
                heater_fluid_bulk_temp).unwrap();

            let generic_advection_interaction = 
            HeatTransferInteractionType::new_advection_interaction(
                mass_flowrate,
                heater_therminol_avg_density,
                heater_therminol_avg_density,
            );
            {
                // prints therminol temperature 
                let heater_surf_temp_degc: Vec<f64> = heater_surface_array_temp
                    .iter().map(
                        |&temperature|{
                            temperature.get::<degree_celsius>()
                        }
                    ).collect();

                // print surface temperature 
                dbg!(heater_surf_temp_degc);

                let bt_12_temperature: ThermodynamicTemperature = 
                static_mixer_pipe_therminol_clone.get_temperature_vector().unwrap() 
                    .into_iter().last().unwrap();
                // print outlet temperature 
                dbg!(heater_top_head_exit_temperature
                .into_format_args(degree_celsius,uom::fmt::DisplayStyle::Abbreviation));

                // bt_12_temperature, which is actually the output temperature of static 
                // mixer 10
                dbg!(bt_12_temperature
                .into_format_args(degree_celsius,uom::fmt::DisplayStyle::Abbreviation));

                //// print therminol temperature 
                //dbg!("Therminol Array Temp: ", therminol_array_temperature);

                //// print twisted tape temperature 
                //dbg!("twisted tape Temp: 
                //note: conduction occurs, so first node is hotter\n 
                //than the therminol fluid", twisted_tape_temperature);

                // print simulation time
                // dbg diagnostics probably not the cause of mem leaks
                dbg!(simulation_time);
            }

            // make axial connections to BCs 
            //
            // note: need to speed up this part, too slow

            heater_bottom_head_bare.therminol_array.link_to_back(
                &mut inlet_bc,
                generic_advection_interaction
            ).unwrap();

            heater_v2_bare.therminol_array.link_to_back(
                &mut heater_bottom_head_bare.therminol_array,
                generic_advection_interaction
            ).unwrap();

            heater_v2_bare.therminol_array.link_to_front(
                &mut heater_top_head_bare.therminol_array,
                generic_advection_interaction
            ).unwrap();


            heater_top_head_bare.therminol_array.link_to_front(
                &mut static_mixer_mx_10_object.therminol_array,
                generic_advection_interaction
            ).unwrap();

            static_mixer_mx_10_object.therminol_array.link_to_front(
                &mut static_mixer_mx_10_pipe.therminol_array,
                generic_advection_interaction
            ).unwrap();

            static_mixer_mx_10_pipe.therminol_array.link_to_front(
                &mut outlet_bc,
                generic_advection_interaction
            ).unwrap();


            // make other connections by spawning a new thread 
            // this is the parallel version
            let heater_2_join_handle: JoinHandle<HeaterVersion2Bare> 
            = heater_v2_bare.
                lateral_connection_thread_spawn(
                    mass_flowrate,
                    heater_power);

            let heater_bottom_join_handle: JoinHandle<HeaterTopBottomHead> 
            = heater_bottom_head_bare. 
                lateral_connection_thread_spawn(
                    mass_flowrate);

            let heater_top_head_join_handle = 
            heater_top_head_bare.lateral_connection_thread_spawn(
                mass_flowrate);


            let static_mixer_join_handle = 
            static_mixer_mx_10_object.lateral_connection_thread_spawn(
                mass_flowrate);

            let static_mixer_pipe_join_handle = 
            static_mixer_mx_10_pipe.lateral_connection_thread_spawn(
                mass_flowrate);

            // link struct supports to ambient air
            // axially 

            structural_support_heater_bottom_head. 
                support_array.link_to_front(
                    &mut ambient_air_temp_bc,
                    support_conductance_interaction
                ).unwrap();

            structural_support_heater_top_head. 
                support_array.link_to_front(
                    &mut ambient_air_temp_bc,
                    support_conductance_interaction
                ).unwrap();

            structural_support_mx_10.support_array.link_to_front(
                &mut ambient_air_temp_bc,
                support_conductance_interaction
            ).unwrap();


            static_mixer_mx_10_object = static_mixer_join_handle.join().unwrap();
            static_mixer_mx_10_pipe = static_mixer_pipe_join_handle.join().unwrap();
            heater_v2_bare = heater_2_join_handle.join().unwrap();
            heater_bottom_head_bare = heater_bottom_join_handle.join().unwrap();
            heater_top_head_bare = heater_top_head_join_handle.join().unwrap();


            // link struct supports to heater top/bottom heads
            structural_support_heater_top_head.
                support_array.link_to_back(
                    &mut heater_top_head_bare.steel_shell,
                    support_conductance_interaction
                ).unwrap();
            structural_support_heater_bottom_head. 
                support_array.link_to_back(
                    &mut heater_bottom_head_bare.steel_shell,
                    support_conductance_interaction
                ).unwrap();

            structural_support_mx_10.support_array.link_to_back(
                &mut static_mixer_mx_10_pipe.steel_shell,
                support_conductance_interaction
            ).unwrap();

            // note, the heater top and bottom head area changed 
            // during course of this interaction, so should be okay


            // i will also connect heater shell to the structural support 
            // via the head as in ciet 

            heater_v2_bare.steel_shell.link_to_back(
                &mut heater_bottom_head_bare.steel_shell,
                support_conductance_interaction
            ).unwrap();

            heater_v2_bare.steel_shell.link_to_front(
                &mut heater_top_head_bare.steel_shell,
                support_conductance_interaction
            ).unwrap();

            // probably edit this to include twisted tape conductance
            heater_v2_bare.twisted_tape_interior.link_to_back(
                &mut heater_bottom_head_bare.twisted_tape_interior,
                support_conductance_interaction
            ).unwrap();

            heater_v2_bare.twisted_tape_interior.link_to_front(
                &mut heater_top_head_bare.twisted_tape_interior,
                support_conductance_interaction
            ).unwrap();

            // now link it laterally to ambient temperatures
            let struct_support_top_head_join_handle = 
            structural_support_heater_top_head.lateral_connection_thread_spawn();
            let structural_support_heater_bottom_head_join_handle = 
            structural_support_heater_bottom_head.lateral_connection_thread_spawn();

            structural_support_mx_10.
                lateral_and_miscellaneous_connections();

            structural_support_heater_top_head = 
                struct_support_top_head_join_handle.join().unwrap();
            structural_support_heater_bottom_head = 
                structural_support_heater_bottom_head_join_handle.join().unwrap();

            //// calculate timestep (serial method)
            //heater_v2_bare.advance_timestep(
            //    timestep);

            // calculate timestep (thread spawn method, parallel) 

            let heater_2_join_handle: JoinHandle<HeaterVersion2Bare> 
            = heater_v2_bare.advance_timestep_thread_spawn(
                timestep);

            let heater_bottom_join_handle: JoinHandle<HeaterTopBottomHead> 
            = heater_bottom_head_bare. 
                advance_timestep_thread_spawn(
                    timestep);

            let heater_top_head_join_handle = 
            heater_top_head_bare.advance_timestep_thread_spawn(
                timestep);

            let static_mixer_join_handle = 
            static_mixer_mx_10_object.advance_timestep_thread_spawn(
                timestep);

            let static_mixer_pipe_join_handle = 
            static_mixer_mx_10_pipe.advance_timestep_thread_spawn(
                timestep);


            let structural_support_heater_bottom_head_join_handle = 
            structural_support_heater_bottom_head.
                advance_timestep_thread_spawn(timestep);
            let structural_support_heater_top_head_join_handle = 
            structural_support_heater_top_head.
                advance_timestep_thread_spawn(timestep);

            structural_support_mx_10._advance_timestep(
                timestep);

            structural_support_heater_bottom_head 
                =  structural_support_heater_bottom_head_join_handle.join().unwrap();
            structural_support_heater_top_head 
                =  structural_support_heater_top_head_join_handle.join().unwrap();


            static_mixer_mx_10_object = static_mixer_join_handle.join().unwrap();
            static_mixer_mx_10_pipe = static_mixer_pipe_join_handle.join().unwrap();
            heater_v2_bare = heater_2_join_handle.join().unwrap();
            heater_bottom_head_bare = heater_bottom_join_handle.join().unwrap();
            heater_top_head_bare = heater_top_head_join_handle.join().unwrap();


            simulation_time += timestep;

            let time_taken_for_calculation_loop = loop_time.elapsed().unwrap()
            - loop_time_start;

            dbg!(time_taken_for_calculation_loop);

        }

    });

    main_loop.join().unwrap();

    let thread_sleep = false;

    if thread_sleep {
        let ten_seconds = time::Duration::from_millis(10000);

        thread::sleep(ten_seconds);
    }


    // once simulation completed, write data


    //todo!("haven't coded csv writing file")

}


