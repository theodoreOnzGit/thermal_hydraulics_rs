
use std::f64::consts::PI;
use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::data_enum_structs::DataUserSpecifiedConvectionResistance;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::{link_heat_transfer_entity, HeatTransferInteractionType};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::SolidMaterial::SteelSS304L;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::{Material, LiquidMaterial};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, OuterDiameterThermalConduction, SurfaceArea, SingleCVNode, CVType, BCType};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::CVType::SingleCV;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::density::density;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy};



use uom::si::angle::radian;
use uom::si::angular_velocity::radian_per_second;
use uom::si::area::square_centimeter;
use uom::si::f64::*;
use uom::si::length::centimeter;
use uom::si::power::{watt, kilowatt};
use uom::si::temperature_interval::degree_celsius as interval_deg_c;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

/// In this test, we have a one dimensional representation of the 
/// CIET heater v1.0 
///
/// Electrical Heat
/// -------------> (solid shell) 
///                     | 
///                     | 
///                     |  (thermal resistance)
///                     | 
///
/// --------> (well mixed fluid volume) ----------->
///                 T_fluid
/// Fluid in                                Fluid out 
/// T_in                                    T_fluid
///
///
/// The outlet temperature is same as the fluid inlet temperature
///
/// I intend to generate a Bode Plot from this using frequency 
/// response data
///
/// I can use either a PRBS signal or MFBS (multifrequency binary signal)
/// to perturb this system. For now, it seems MFBS is an easier choice 
/// as PRBS signal generators in rust are a little harder to do 
///
/// So for MFBS, we can superimpose about 10 sine waves, each with 
/// its own frequency, sample the time to obtain a value 
/// If the value is greater then zero, then the output is 1 
/// otherwise, the output is zero
///
/// Use that to determine temperature input signals over the timestep
/// omega is 0.002301 rad/s to 23.01 rad/s sampled at five frequencies 
/// in here
/// https://fhr.nuc.berkeley.edu/wp-content/uploads/2017/11/TH-Report-UCBTH-17-002-Full-Nov-Update.pdf
/// 
/// equivalently, this is 0.000366 Hz to 3.66 Hz oscillation frequency
/// in matlab, the Bode plot data is from 0.001 rad/s to 10 rad/s
/// 
/// perhaps I shall just use these five angular frequencies to sample 
/// points, of course, using MFBS allows me to sample other signals 
/// using a weaker signal power too
///
///
///
#[test]
pub fn one_dimension_ciet_heater_v_1_0_test(){

    // okay, let's make two control volumes 
    // one cylinder and then the other a shell
    //
    // cylinder needs diameter and z 
    // shell needs id, od and z
    
    let therminol_cylinder: HeatTransferEntity;
    let steel_shell: HeatTransferEntity;

    // need two boundary conditions 

    let inlet_const_temp = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ThermodynamicTemperature::new::<degree_celsius>(80.0)
        ));

    let outlet_zero_heat_flux = HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatAddition(Power::new::<watt>(0.0))
    );

    // the two types of HeatTransferInteractionType are 
    // advection and convection resistance
    //
    // 2007 square_centimeter
    // and 607 watt_per_square_meter_kelvin

    let convection_data = DataUserSpecifiedConvectionResistance { 
        surf_area: SurfaceArea::from(
            Area::new::<square_centimeter>(2007_f64)
        ),
        heat_transfer_coeff: 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(607_f64),
    };

    let convection_resistance = HeatTransferInteractionType::
        UserSpecifiedConvectionResistance(
            convection_data
        );



    // timestep settings


    let max_time: Time = Time::new::<second>(20.0);
    let max_time_ptr = Arc::new(max_time);

    // this is the calculation loop
    let calculation_loop = move || {
        let mut current_time_simulation_time = Time::new::<second>(0.0);

        let max_time_ptr_in_loop = max_time_ptr;
        // we are sampling at about 10 Hz
        // so the nyquist frequency is about 5 Hz 
        // this is because the highest frequency is about 3.66 Hz
        let mut timestep_value = Time::new::<second>(0.1);
        
        while current_time_simulation_time <= *max_time_ptr_in_loop {

            // todo calculation steps

            // advection bc, so at boundary condition, therminol flows in at 
            // 0.18 kg/s at 80 c 

            let therminol = Material::Liquid(LiquidMaterial::TherminolVP1);

            let therminol_inlet_temperature = 
            ThermodynamicTemperature::new::<degree_celsius>(80.0);
            let atmospheric_pressure = 
            Pressure::new::<atmosphere>(1.0);

            let therminol_inlet_density = density(
                therminol,
                therminol_inlet_temperature,
                atmospheric_pressure
            ).unwrap();

            // i can also calculate the densities of each cv



            // todo csv data writing

            // add the timestep
            current_time_simulation_time += timestep_value;
        }
    };

    let calculation_thread = thread::spawn(calculation_loop);
    
    calculation_thread.join().unwrap();

    // done!
    return ();
}

fn mfbs_poresky_2017_power_signal(simulation_time: Time) -> Power {

    // get simulation time in seconds 

    // obtain power, which is 9 kW plus 1kW amplitude 
    
    let steady_state_power = Power::new::<kilowatt>(9.0);
    let power_amplitude = Power::new::<kilowatt>(1.0);

    let mut power_vector: Vec<Power> = vec![];

    let mut angular_frequency_vector: Vec<AngularVelocity>
    = vec![];

    angular_frequency_vector.push(
        AngularVelocity::new::<radian_per_second>(0.002301));
    angular_frequency_vector.push(
        AngularVelocity::new::<radian_per_second>(0.02301));
    angular_frequency_vector.push(
        AngularVelocity::new::<radian_per_second>(0.2301));
    angular_frequency_vector.push(
        AngularVelocity::new::<radian_per_second>(2.301));
    angular_frequency_vector.push(
        AngularVelocity::new::<radian_per_second>(23.01));

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
        power_analog_signal += *power_ptr;
    }

    // if power signal is less than 9kW, then return 8kW otherwise 
    // 10 kW
    //
    // kinda bloated code, but it does have the analog signal 
    // if i want

    if power_analog_signal.get::<kilowatt>() < 9.0 {
        steady_state_power - power_amplitude
    } else {
        steady_state_power + power_amplitude
    }


}

