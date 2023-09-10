use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::SystemTime;

use csv::Writer;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, RadialCylindricalThicknessThermalConduction, InnerDiameterThermalConduction};
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::one_dimension_fluid_array::FluidArray;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::one_dimension_solid_array::SolidColumn;
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::{SolidMaterial, Material};
use thermal_hydraulics_rs::heat_transfer_lib::
thermophysical_properties::LiquidMaterial;



use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::dynamic_viscosity::dynamic_viscosity;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::prandtl::liquid_prandtl;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use uom::si::angle::radian;
use uom::si::angular_velocity::radian_per_second;
use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::{centimeter, meter};
use uom::si::mass_rate::kilogram_per_second;
use uom::si::power::{watt, kilowatt};
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;




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

    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::
    heat_transfer_interactions::data_enum_structs::DataAdvection;
    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::
    heat_transfer_interactions::data_enum_structs::DataUserSpecifiedConvectionResistance;
    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::
    heat_transfer_interactions::HeatTransferInteractionType;
    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::
    heat_transfer_interactions::link_heat_transfer_entity;
    use thermal_hydraulics_rs::heat_transfer_lib::
    thermophysical_properties::Material;

    use uom::si::length::inch;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
    ::heat_transfer_entities::BCType::*;
    use uom::si::thermodynamic_temperature::kelvin;

    use ndarray::*;
    // okay, let's make two control volumes 
    // one cylinder and then the other a shell
    //
    // cylinder needs diameter and z 
    // shell needs id, od and z
    let therminol = Material::Liquid(LiquidMaterial::TherminolVP1);
    let steel = Material::Solid(SolidMaterial::SteelSS304L);
    let id = Length::new::<meter>(0.0381);
    let od = Length::new::<meter>(0.04);
    let _inner_tube_od = Length::new::<centimeter>(3.175);
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
    let node_length: Length = heated_length/number_of_nodes as f64;

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

    let therminol_array: HeatTransferEntity = 
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
    ).into();

    // now the outer steel array

    let steel_shell_array: HeatTransferEntity = 
    SolidColumn::new_cylindrical_shell(
        heated_length,
        id,
        od,
        initial_temperature,
        atmospheric_pressure,
        SolidMaterial::SteelSS304L,
        6
    ).into();

    // move arrays to Arc ptrs 

    let fluid_array_ptr = Arc::new(Mutex::new(
        therminol_array
    ));

    let solid_array_ptr = Arc::new(Mutex::new( 
        steel_shell_array 
    ));

    // clone ptrs for moving into loop 
    
    let fluid_array_ptr_for_loop = fluid_array_ptr.clone();
    let solid_array_ptr_for_loop = solid_array_ptr.clone();

    // time settings
    let max_time: Time = Time::new::<second>(0.02);
    let max_time_ptr = Arc::new(max_time);

    let calculation_time_elapsed = SystemTime::now();
    let calculation_loop = move || {
        // timestep settings

        let mut current_time_simulation_time = Time::new::<second>(0.0);
        let timestep = Time::new::<second>(0.01);


        //csv writer

        let mut time_wtr = Writer::from_path("fluid_node_calc_time_profile.csv")
            .unwrap();

        time_wtr.write_record(&["loop_calculation_time_nanoseconds",
            "mutex_lock_time_ns",
            "node_connection_time_ns",
            "data_record_time_ns",
            "timestep_advance_time_ns",])
            .unwrap();

        let mut temp_profile_wtr = Writer::from_path(
            "fluid_node_temp_profile.csv")
            .unwrap();

        // this is code for writing the array of required temperatures
        {

            // I want the mid node length of this temperature


            let half_node_length: Length = 0.5 * node_length;

            let mut header_vec: Vec<String> = vec![];

            header_vec.push("simulation_time_seconds".to_string());
            header_vec.push("elapsed_time_seconds".to_string());

            for index in 0..number_of_nodes {

                let mid_node_length: Length = 
                index as f64 * node_length + half_node_length;

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

        }

        // main calculation loop

        while current_time_simulation_time <= *max_time_ptr.deref(){

            let arc_mutex_lock_start = SystemTime::now();
            // obtain arc mutex locks 


            let mut therminol_entity_ptr_in_loop = 
            fluid_array_ptr_for_loop.lock().unwrap();

            let mut steel_entity_ptr_in_loop 
            = solid_array_ptr_for_loop.lock().unwrap();
            
            let arc_mutex_lock_elapsed_ms = 
            arc_mutex_lock_start.elapsed().unwrap().as_millis();


            // next, obtain average temperature 
            // average by volume for both fluid and solid vec
            // so we can get an average conductance value
            //
            let node_connection_start = SystemTime::now();

            // now I'm going to clone the heat transfer entity 
            // as a fluid array and solid column respectively to 
            // get the underlying bulk temp

            let mut therminol_array_clone: FluidArray = 
            therminol_entity_ptr_in_loop.deref_mut().clone().try_into().unwrap();
            let mut steel_array_clone: SolidColumn = 
            steel_entity_ptr_in_loop.deref_mut().clone().try_into().unwrap();

            let fluid_average_temp: ThermodynamicTemperature = 
            therminol_array_clone.get_bulk_temperature().unwrap();

            let solid_average_temp: ThermodynamicTemperature = 
            steel_array_clone.get_bulk_temperature().unwrap();

            // given these two, we can calculate an average conductance 
            // value across all solid-fluid boundaries. 
            //
            // This is a time saving measure, rather than calculating 
            // all the conductances node by node

            // now, we make a cylindrical conductance interaction with 
            // fluid inside

            // we'll need the thickness 

            let radial_thickness_thermal_conduction = 0.5*(od - id);
            let radial_thickness_thermal_conduction: 
            RadialCylindricalThicknessThermalConduction = 
            radial_thickness_thermal_conduction.into();

            let inner_diameter_thermal_conduction: InnerDiameterThermalConduction 
            = id.into();


            // need to change this the nusselt correlation for the actual 
            // test
            // heater 2.0
            let h_to_therminol: HeatTransfer = 
            _heat_transfer_coefficient_ciet_v_2_0(
                therminol_mass_flowrate,
                fluid_average_temp,
                atmospheric_pressure);

            let therminol_steel_conductance_interaction: HeatTransferInteractionType
            = HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(
                    (steel, 
                    radial_thickness_thermal_conduction,
                    solid_average_temp,
                    atmospheric_pressure),
                    (h_to_therminol,
                    inner_diameter_thermal_conduction,
                    node_length.clone().into())
            );

            // now based on conductance interaction, 
            // we can obtain thermal conductance, the temperatures 
            // and pressures don't really matter
            //
            // this is because all the thermal conductance data 
            // has already been loaded into the thermal conductance 
            // interaction object

            let therminol_steel_nodal_thermal_conductance: ThermalConductance = 
            therminol_steel_conductance_interaction.try_get_thermal_conductance(
                fluid_average_temp,
                solid_average_temp,
                atmospheric_pressure,
                atmospheric_pressure,
            ).unwrap();

            // now, create conductance vector  for heater v2.0

            let mut steel_therminol_conductance_vector: Array1<ThermalConductance> = 
            Array::zeros(number_of_nodes);



            steel_therminol_conductance_vector.fill(therminol_steel_nodal_thermal_conductance);
            // advance timestep
            current_time_simulation_time += timestep;

            
        }
    };

    let main_thread_handle = thread::spawn(calculation_loop);

    main_thread_handle.join().unwrap();


    return ();

}



fn _mfbs_power_signal_logspace_custom(simulation_time: Time) -> Power {

    // get simulation time in seconds 

    // obtain power, which is 9 kW plus 1kW amplitude 
    
    let steady_state_power = Power::new::<kilowatt>(9.0);
    let power_amplitude = Power::new::<kilowatt>(1.0);

    let mut power_vector: Vec<Power> = vec![];

    let angular_frequency_vector_values: Vec<f64> 
    = _logspace(0.001, 1_f64, 15);

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
fn _logspace(start: f64, end: f64, n: usize) -> Vec<f64> {
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
fn _ciet_heater_v_2_0_nusselt_number(reynolds:Ratio, 
    prandtl:Ratio) -> Ratio {

    let reynolds_power_0_836 = reynolds.value.powf(0.836);
    let prandtl_power_0_333 = prandtl.value.powf(0.333333333333333);

    Ratio::new::<ratio>(
    0.04179 * reynolds_power_0_836 * prandtl_power_0_333)

}
#[inline]
fn _ciet_heater_v_2_0_reynolds_nunber(mass_flowrate: MassRate,
    mu: DynamicViscosity) -> Ratio {

    // Re = m* D_H/ A_{XS}/mu
    let hydraulic_diameter = Length::new::<meter>(0.01467);
    let flow_area = Area::new::<square_meter>(0.00105);

    mass_flowrate*hydraulic_diameter/mu/flow_area
}

fn _heat_transfer_coefficient_ciet_v_2_0(mass_flowrate: MassRate,
    therminol_temperature: ThermodynamicTemperature,
    pressure: Pressure) -> HeatTransfer {

    // let's calculate mu and k 

    let therminol = Material::Liquid(LiquidMaterial::TherminolVP1);
    let mu: DynamicViscosity = dynamic_viscosity(therminol,
        therminol_temperature,
        pressure).unwrap();

    let k: ThermalConductivity = thermal_conductivity(
        therminol,
        therminol_temperature,
        pressure).unwrap();


    let reynolds: Ratio = _ciet_heater_v_2_0_reynolds_nunber(
        mass_flowrate, mu);

    let prandtl: Ratio = liquid_prandtl(
        therminol,
        therminol_temperature,
        pressure).unwrap();

    let nusselt: Ratio = _ciet_heater_v_2_0_nusselt_number(
        reynolds,
        prandtl);

    let hydraulic_diameter = Length::new::<meter>(0.01467);

    let heat_transfer_coeff: HeatTransfer = 
    nusselt * k / hydraulic_diameter;

    heat_transfer_coeff

}
