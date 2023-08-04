use std::f64::consts::PI;
use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;

use peroxide::prelude::erfc;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::data_enum_structs::DataUserSpecifiedConvectionResistance;
use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::{link_heat_transfer_entity, HeatTransferInteractionType, calculate_timescales_for_heat_transfer_entity};
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

    // let's do the thread spawn for the calculated solution before 
    // the analytical solution

    // before we start, we need the copper thermal_diffusivity

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

    let copper_initial_enthalpy = specific_enthalpy(
        copper, 
        copper_initial_temperature, 
        pressure)?;

    let node_length: Length = Length::new::<centimeter>(4.0);
    // let's make the up to 0.40 m of nodes, 
    // we have 10 nodes in all, so each node has
    // about 4cm of thermal resistance between the nodes
    //
    // then specify each node is made of copper, 
    // it shall be a 1D kind of conduction node
    // for 1D conduction type nodes, we just take 1m^2 area as a basis 
    // and apply adiabatic BC to the surface normal 
    // 
    // so I'll need to make 1 BC and 10 control volumes
    // simulate the transient temperatures for up to 20s

    let single_cv_node_1 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_2 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_3 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_4 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_5 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_6 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_7 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_8 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_9 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;
    let single_cv_node_10 = SingleCVNode::new_one_dimension_volume(
        node_length, copper, copper_initial_temperature, 
        pressure)?;

    let single_cv_1_ptr = Arc::new(
        Mutex::new(single_cv_node_1)
    );
    let single_cv_2_ptr = Arc::new(
        Mutex::new(single_cv_node_2)
    );
    let single_cv_3_ptr = Arc::new(
        Mutex::new(single_cv_node_3)
    );
    let single_cv_4_ptr = Arc::new(
        Mutex::new(single_cv_node_4)
    );
    let single_cv_5_ptr = Arc::new(
        Mutex::new(single_cv_node_5)
    );
    let single_cv_6_ptr = Arc::new(
        Mutex::new(single_cv_node_6)
    );
    let single_cv_7_ptr = Arc::new(
        Mutex::new(single_cv_node_7)
    );
    let single_cv_8_ptr = Arc::new(
        Mutex::new(single_cv_node_8)
    );
    let single_cv_9_ptr = Arc::new(
        Mutex::new(single_cv_node_9)
    );
    let single_cv_10_ptr = Arc::new(
        Mutex::new(single_cv_node_10)
    );

    let copper_surface_temperature_boundary_condition = 
    HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(boundary_condition_temperature));

    let surf_temp_ptr = Arc::new(
        Mutex::new(copper_surface_temperature_boundary_condition)
    );

    // time settings
    let max_time: Time = Time::new::<second>(20.0);

    let timestep: Time = Time::new::<second>(0.1);
    let timestep_ptr = Arc::new(
        Mutex::new(timestep)
    );
    let max_time_ptr = Arc::new(max_time);

    let calculation_loop = move || {

        let mut single_cv_1_in_loop = single_cv_1_ptr.lock().unwrap();
        let mut single_cv_2_in_loop = single_cv_2_ptr.lock().unwrap();
        let mut single_cv_3_in_loop = single_cv_3_ptr.lock().unwrap();
        let mut single_cv_4_in_loop = single_cv_4_ptr.lock().unwrap();
        let mut single_cv_5_in_loop = single_cv_5_ptr.lock().unwrap();
        let mut single_cv_6_in_loop = single_cv_6_ptr.lock().unwrap();
        let mut single_cv_7_in_loop = single_cv_7_ptr.lock().unwrap();
        let mut single_cv_8_in_loop = single_cv_8_ptr.lock().unwrap();
        let mut single_cv_9_in_loop = single_cv_9_ptr.lock().unwrap();
        let mut single_cv_10_in_loop = single_cv_10_ptr.lock().unwrap();

        let mut surf_temp_bc_in_loop = surf_temp_ptr.lock().unwrap();

        let mut timestep_in_loop = timestep_ptr.lock().unwrap();
        let max_time_ptr_in_loop = max_time_ptr;

        use csv::Writer;
        let mut wtr = Writer::from_path("semi_infinite_simulated_values.csv")
            .unwrap();

        wtr.write_record(&["time_seconds",
            "0cm_temperautre_kelvin",
            "2cm_temperature_kelvin",
            "6cm_temperature_kelvin",
            "10cm_temperature_kelvin",
            "14cm_temperature_kelvin",
            "18cm_temperature_kelvin",
            "22cm_temperature_kelvin",
            "26cm_temperature_kelvin",
            "30cm_temperature_kelvin",
            "34cm_temperature_kelvin",
            "38cm_temperature_kelvin",
        ]).unwrap();

        // let's establish interactions between each of the nodes
        //
        // the first node would have a 2cm thermal resistance as it 
        // is closest to the wall 
        //
        //
        // | 
        // | 
        // | 2cm                4cm                 4cm
        // -------- * ------------------------- * ---------------
        // |        node 1                      node 2 
        // | 
        // | 
        // 

        let first_node_thermal_resistance = 
        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(copper,
                Length::new::<centimeter>(2.0).into());


        let subsequent_node_thermal_resistance = 
        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(copper,
                Length::new::<centimeter>(4.0).into());


        let mut current_time_simulation_time = Time::new::<second>(0.0);


        while current_time_simulation_time <= *max_time_ptr_in_loop {

            // first let's link the heat transfer entities 


            // first node is very important, we have BC and CV linkage
            link_heat_transfer_entity(&mut surf_temp_bc_in_loop,
                &mut single_cv_1_in_loop,
                first_node_thermal_resistance).unwrap();

            // subsequent nodes have similar linkages
            link_heat_transfer_entity(&mut single_cv_1_in_loop,
                &mut single_cv_2_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_2_in_loop,
                &mut single_cv_3_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_3_in_loop,
                &mut single_cv_4_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_4_in_loop,
                &mut single_cv_5_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_5_in_loop,
                &mut single_cv_6_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_6_in_loop,
                &mut single_cv_7_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_7_in_loop,
                &mut single_cv_8_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_8_in_loop,
                &mut single_cv_9_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            link_heat_transfer_entity(&mut single_cv_9_in_loop,
                &mut single_cv_10_in_loop,
                subsequent_node_thermal_resistance).unwrap();

            // now let's capture the temperature data first 

            let bc_temperature: ThermodynamicTemperature 
            = boundary_condition_temperature;

            let cv_1_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_1_in_loop.deref_mut()).unwrap();

            let cv_2_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_2_in_loop.deref_mut()).unwrap();

            let cv_3_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_3_in_loop.deref_mut()).unwrap();

            let cv_4_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_4_in_loop.deref_mut()).unwrap();

            let cv_5_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_5_in_loop.deref_mut()).unwrap();

            let cv_6_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_6_in_loop.deref_mut()).unwrap();

            let cv_7_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_7_in_loop.deref_mut()).unwrap();

            let cv_8_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_8_in_loop.deref_mut()).unwrap();

            let cv_9_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_9_in_loop.deref_mut()).unwrap();

            let cv_10_temperature = 
            HeatTransferEntity::temperature( 
                single_cv_10_in_loop.deref_mut()).unwrap();

            let time_string = current_time_simulation_time.value.to_string();

            wtr.write_record(&[time_string,
                bc_temperature.value.to_string(),
                cv_1_temperature.value.to_string(),
                cv_2_temperature.value.to_string(),
                cv_3_temperature.value.to_string(),
                cv_4_temperature.value.to_string(),
                cv_5_temperature.value.to_string(),
                cv_6_temperature.value.to_string(),
                cv_7_temperature.value.to_string(),
                cv_8_temperature.value.to_string(),
                cv_9_temperature.value.to_string(),
                cv_10_temperature.value.to_string(), ])
                .unwrap();

            // now we need to update the timestep 
            // we'll just use the cv-bc timestep because that has 
            // the smallest lengthscale, should be the shortest

            let timestep_from_api = 
            calculate_timescales_for_heat_transfer_entity(
                &mut surf_temp_bc_in_loop,
                &mut single_cv_1_in_loop,
                first_node_thermal_resistance).unwrap();


            let timestep_value = timestep_from_api;
            // update timestep value

            *timestep_in_loop.deref_mut() = timestep_value;

            HeatTransferEntity::advance_timestep(
                single_cv_1_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();

            HeatTransferEntity::advance_timestep(
                single_cv_2_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_3_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_4_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_5_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_6_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_7_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_8_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_9_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();
            HeatTransferEntity::advance_timestep(
                single_cv_10_in_loop.deref_mut(),
                *timestep_in_loop).unwrap();

            current_time_simulation_time += *timestep_in_loop.deref();
        }
        wtr.flush().unwrap();
    };

    let calculation_thread = thread::spawn(calculation_loop);



    // then let's do the analytical solution
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

    calculation_thread.join().unwrap();
    // this setup is meant to be emulated using control volumes with 
    // some thermal resistances between them

    //todo!("need to do 4 control vol and 1 BC for 1d transient conduction");
    Ok(())
}
