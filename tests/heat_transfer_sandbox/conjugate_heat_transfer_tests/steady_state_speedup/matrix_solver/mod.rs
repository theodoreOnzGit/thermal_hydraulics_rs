use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::SystemTime;

use csv::Writer;
use ndarray::Array1;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::one_dimension_fluid_array::FluidArray;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::data_enum_structs::{DataUserSpecifiedConvectionResistance, DataAdvection};
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::{link_heat_transfer_entity, HeatTransferInteractionType};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::SolidMaterial::{self};
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::dynamic_viscosity::dynamic_viscosity;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::prandtl::liquid_prandtl;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::{Material, LiquidMaterial};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, 
    SingleCVNode, CVType, BCType, InnerDiameterThermalConduction, OuterDiameterThermalConduction, RadialCylindricalThicknessThermalConduction, CylinderLengthThermalConduction};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::CVType::SingleCV;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::density::density;



use uom::si::angle::radian;
use uom::si::angular_velocity::radian_per_second;
use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::{centimeter, meter, inch};
use uom::si::mass_rate::kilogram_per_second;
use uom::si::power::{watt, kilowatt};
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

use ndarray::*;



/// In this test, we have a nodalised representation of the 
/// this is 8 nodes in the axial direction and 2 nodes for metal 
/// in the radial direction
///
/// This is for heater v2.0
///
/// I realised that we needed a major speedup because in debug mode 
/// a 2.3ms timestep was taking ~ 602 ms to calculate
///
/// This was unacceptable, especially when we need about 30 plus 
/// components for in forced convection loop in CTAH plus heater branch
///
/// Hence there is a need not only to speed up calculations, but also 
/// to quickly initialise nodalised control volumes. For example, 
/// the heater has 8 nodes
///
/// For this test, we shall design a heater with 8 nodes axially 
/// one fluid node 
/// and one steel node
/// 
/// (ambient air)
///
/// | 
/// | 
/// | 
/// 
/// (steel shell) ---- (heater power)
///
/// | 
/// |  
/// | 
///
/// (fluid node)  --- (subsequent nodes, advection)
///
/// each of these nodes will be represented by matrices
///
/// 
///
#[test]
//#[ignore = "data collected"]
pub fn matrix_calculation_initial_test(){


    // okay, let's make two control volumes 
    // one cylinder and then the other a shell
    //
    // cylinder needs diameter and z 
    // shell needs id, od and z
    let therminol = Material::Liquid(LiquidMaterial::TherminolVP1);
    let steel = Material::Solid(SolidMaterial::SteelSS304L);
    let id = Length::new::<meter>(0.0381);
    let od = Length::new::<meter>(0.04);
    let inner_tube_od = Length::new::<centimeter>(3.175);
    // z is heated length
    let _total_length = Length::new::<meter>(1.983333);
    let heated_length = Length::new::<meter>(1.676);
    let initial_temperature = ThermodynamicTemperature::new::
        <degree_celsius>(80.0);
    let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

    let flow_area = Area::new::<square_meter>(0.00105);
    let number_of_nodes: usize = 8;
    let ambient_air_temp = ThermodynamicTemperature::new::<
        degree_celsius>(21.67);

    let heater_steady_state_power = Power::new::<kilowatt>(8.0);
    let therminol_mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);

    // the form loss of the pipe is meant to help construct the inner 
    // pipe object, that is to help it return the nusselt number 
    // that is expected for this pipe
    // 
    // However, one is free to specify whatever nusselt number works
    let dummy_pipe_form_loss = Ratio::new::<ratio>(0.1);

    // heater is inclined 90 degrees upwards, not that this is 
    // particularly important for this scenario
    let pipe_incline_angle = Angle::new::<uom::si::angle::degree>(90.0);
    // now, to construct the heater
    //
    // I'll construct the inner fluid tube first 

    let therminol_array: FluidArray = 
    FluidArray::new_odd_shaped_pipe(
        heated_length,
        flow_area,
        initial_temperature,
        atmospheric_pressure,
        SolidMaterial::SteelSS304L,
        LiquidMaterial::TherminolVP1,
        dummy_pipe_form_loss,
        6,
        pipe_incline_angle
    );




    // move the fluid temp arrays into arc ptrs with mutex lock 
    // as well as the single cvs and such 


    // time settings
    let max_time: Time = Time::new::<second>(0.02);
    let max_time_ptr = Arc::new(max_time);

    let calculation_time_elapsed = SystemTime::now();
    let calculation_loop = move || {
        // derference ptrs 

        let mut current_time_simulation_time = Time::new::<second>(0.0);
        let timestep = Time::new::<second>(0.01);

        while current_time_simulation_time <= *max_time_ptr.deref(){

            let arc_mutex_lock_start = SystemTime::now();
            // obtain arc mutex locks 


            // advance timestep
            current_time_simulation_time += timestep;
            
            // time for arc mutex derefs
            // this may cause race conditions
            //let arc_mutex_lock_end = SystemTime::now();

            //let arc_mutex_lock_elapsed_ms = 
            //arc_mutex_lock_start.elapsed().unwrap().as_millis()
            //- arc_mutex_lock_end.elapsed().unwrap().as_millis();
            
            let arc_mutex_lock_elapsed_ms = 
            arc_mutex_lock_start.elapsed().unwrap().as_millis();
            // next, obtain average temperature 
            // average by volume for both fluid and solid vec
            //



            
        }
    };

    let main_thread_handle = thread::spawn(calculation_loop);

    main_thread_handle.join().unwrap();


    return ();

}



fn mfbs_power_signal_logspace_custom(simulation_time: Time) -> Power {

    // get simulation time in seconds 

    // obtain power, which is 9 kW plus 1kW amplitude 
    
    let steady_state_power = Power::new::<kilowatt>(9.0);
    let power_amplitude = Power::new::<kilowatt>(1.0);

    let mut power_vector: Vec<Power> = vec![];

    let angular_frequency_vector_values: Vec<f64> 
    = logspace(0.001, 1_f64, 15);

    let mut angular_frequency_vector: Vec<AngularVelocity>
    = vec![];

    for angular_freq_val in angular_frequency_vector_values {

        angular_frequency_vector.push(
            AngularVelocity::new::<radian_per_second>(angular_freq_val));
    }


    for angular_frequency_ptr in angular_frequency_vector.iter() {
        let phase: Angle = (*angular_frequency_ptr * simulation_time).into();

        let phase_angle_radian_value: f64 = phase.get::<radian>();

        let phase_angle_sine: f64 = phase_angle_radian_value.sin();

        let power_component = steady_state_power + 
        power_amplitude * phase_angle_sine;

        power_vector.push(power_component);

    }

    // sum all power components 
    let mut power_analog_signal = Power::new::<watt>(0.0);

    for power_ptr in power_vector.iter() {
        // get average
        power_analog_signal += *power_ptr/power_vector.len() as f64;
    }
    
    // change to mfbs
    if power_analog_signal.get::<kilowatt>() < 9.0 {
        steady_state_power - power_amplitude
    } else {
        steady_state_power + power_amplitude
    }

}

/// generated by perplexity AI,
/// I still needed some mild correction
fn logspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    let base = 10.0_f64;
    let start_log = start.log10();
    let end_log = end.log10();
    let step = (end_log - start_log) / (n - 1) as f64;
    (0..n)
        .map(|i| base.powf(start_log + i as f64 * step))
        .collect()
}

// nusselt number correlation 
#[inline]
fn ciet_heater_v_2_0_nusselt_number(reynolds:Ratio, 
    prandtl:Ratio) -> Ratio {

    let reynolds_power_0_836 = reynolds.value.powf(0.836);
    let prandtl_power_0_333 = prandtl.value.powf(0.333333333333333);

    Ratio::new::<ratio>(
    0.04179 * reynolds_power_0_836 * prandtl_power_0_333)

}

