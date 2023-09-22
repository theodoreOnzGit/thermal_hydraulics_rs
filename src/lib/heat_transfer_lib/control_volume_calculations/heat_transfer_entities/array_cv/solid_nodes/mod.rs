/// code for solving solid nodes connected to two 
/// adjacent boundaries: 
///
/// an inner boundary 
/// and an outer boundary
pub (in crate) mod shell_solid_node;

/// code for solving solid nodes connected to one adjacent boundary 
pub (in crate) mod core_solid_nodes;


use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::area::square_meter;
use uom::si::length::meter;

use crate::heat_transfer_lib::thermophysical_properties::LiquidMaterial;
use crate::heat_transfer_lib::thermophysical_properties::prandtl::try_get_prandtl;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::try_get_kappa_thermal_conductivity;
use crate::heat_transfer_lib::thermophysical_properties::dynamic_viscosity::try_get_mu_viscosity;
use crate::heat_transfer_lib::thermophysical_properties::Material;


#[test]
//#[ignore="data collected"]
pub fn heater_v_2_0_nodalised_matrix_solver_test(){

    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::SingleCVNode;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::RadialCylindricalThicknessThermalConduction;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::InnerDiameterThermalConduction;
    use crate::heat_transfer_lib::control_volume_calculations:: 
    heat_transfer_entities::CVType::*;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::HeatTransferInteractionType;
    use std::{sync::{Arc, Mutex}, time::SystemTime, ops::{Deref, DerefMut}, thread, f64::consts::PI};
    use csv::Writer;
    use ndarray::*;
    
    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use uom::si::time::second;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::area::square_meter;
    use uom::si::pressure::atmosphere;
    use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
    use uom::si::length::{meter, centimeter, inch};
    use uom::si::power::{watt, kilowatt};

    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::fluid_nodes::core_fluid_node::_advance_timestep_fluid_node_array_pipe_high_peclet_number;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::OuterDiameterThermalConduction;
    use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::get_thermal_conductance_based_on_interaction;
    use crate::heat_transfer_lib::thermophysical_properties::LiquidMaterial;
    use crate::heat_transfer_lib::thermophysical_properties::volumetric_heat_capacity::try_get_rho_cp;
    use crate::heat_transfer_lib::thermophysical_properties::SolidMaterial;
    use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::try_get_h;
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

    let inlet_temperature = ThermodynamicTemperature::new::<degree_celsius>
        (79.12);
    // now, to construct the heater
    // solid fluid interactions, we use a row of 
    // temperature vectors, a front CV and a back CV to 
    // The front and back CV should be identical 
    //
    // the 
    let node_length: Length = heated_length/number_of_nodes as f64;

    let fluid_back_entity: HeatTransferEntity = 
    SingleCVNode::new_odd_shaped_pipe(
        node_length,
        flow_area,
        therminol,
        initial_temperature,
        atmospheric_pressure,
    ).unwrap();

    // get the cv inside, 
    // real ugly match statement, but whatever.

    let fluid_back_cv: SingleCVNode = match fluid_back_entity {
        HeatTransferEntity::ControlVolume(cv_type) => {
            match cv_type {
                SingleCV(cv) => {
                    cv
                },
                ArrayCV(_) => {
                    panic!()
                },
            }
        },
        HeatTransferEntity::BoundaryConditions(_) => {
            panic!()
        },
    };

    let fluid_front_cv: SingleCVNode = fluid_back_cv.clone();
    let total_fluid_volume = flow_area * heated_length;

    // total steel volume is also needed 

    let total_steel_volume = heated_length * 
    (PI * 0.25 * od * od 
    - PI * 0.25 * id * id);


    let mut initial_temperature_array: Array1<ThermodynamicTemperature> = 
    Array::default(number_of_nodes);
    initial_temperature_array.fill(initial_temperature);

    let fluid_temperature_array = initial_temperature_array.clone();

    let mut steel_temperature_array: Array1<ThermodynamicTemperature> = 
    Array::default(number_of_nodes);
    
    steel_temperature_array.fill(initial_temperature);

    let mut ambient_air_temperature_array: Array1<ThermodynamicTemperature> = 
    Array::default(number_of_nodes);

    ambient_air_temperature_array.fill(ambient_air_temp);




    // we'll need solid fluid conductance array, volume fraction array 
    // rho_cp array and q_fraction array 
    //
    // for conductance, what I'd like to do is 
    // take an average temperature for both solid and fluid 
    // and then based on that, obtain a conductance value, thus 
    // reducing the number of times we calculate with conductance
    //
    // however, the solid fluid conductance is periodically updated, 
    // so do it in a while loop, same for rhocp

    let mut q_fraction: Array1<f64> = Array::zeros(number_of_nodes);
    q_fraction.fill(1.0/number_of_nodes as f64);

    let vol_fraction_equal_split: f64 = 1.0/number_of_nodes as f64;
    let mut volume_fraction_array: Array1<f64> = Array::zeros(number_of_nodes);
    volume_fraction_array.fill(vol_fraction_equal_split);

    // move the fluid temp arrays into arc ptrs with mutex lock 
    // as well as the single cvs and such 

    let back_cv_ptr = Arc::new(Mutex::new(
        fluid_back_cv
    ));
    let front_cv_ptr = Arc::new(Mutex::new(
        fluid_front_cv
    ));
    let fluid_temp_at_present_timestep_ptr = Arc::new(Mutex::new(
        fluid_temperature_array
    ));
    let ambient_air_temp_at_present_timestep_ptr = Arc::new(Mutex::new(
        ambient_air_temperature_array
    ));
    let q_fraction_ptr = Arc::new(Mutex::new(
        q_fraction
    ));
    let fluid_vol_fraction_ptr = Arc::new(Mutex::new(
        volume_fraction_array.clone()
    ));
    let solid_vol_fraction_ptr = Arc::new(Mutex::new(
        volume_fraction_array
    ));
    let steel_temperature_array_present_timestep_ptr = Arc::new(Mutex::new(
        steel_temperature_array
    ));

    // clone the ptrs for the loop
    let back_cv_ptr_clone_for_loop = back_cv_ptr.clone();
    let front_cv_ptr_clone_for_loop = front_cv_ptr.clone();
    let fluid_temp_vec_ptr_for_loop = fluid_temp_at_present_timestep_ptr.clone();
    let ambient_air_temp_at_present_timestep_ptr_for_loop = 
    ambient_air_temp_at_present_timestep_ptr.clone();
    let q_fraction_ptr_for_loop = q_fraction_ptr.clone();
    let fluid_vol_fraction_ptr_for_loop = fluid_vol_fraction_ptr.clone();
    let solid_vol_fraction_ptr_for_loop = solid_vol_fraction_ptr.clone();
    let steel_temperature_array_present_timestep_ptr_for_loop = 
    steel_temperature_array_present_timestep_ptr.clone();

    // time settings
    let max_time: Time = Time::new::<second>(200.0);
    let max_time_ptr = Arc::new(max_time);

    let calculation_time_elapsed = SystemTime::now();
    let calculation_loop = move || {
        // derference ptrs 

        let mut current_time_simulation_time = Time::new::<second>(0.0);
        let timestep = Time::new::<second>(0.01);

        // csv writer (heater 2.0)

        let mut wtr = Writer::from_path("matrix_solve_trial_ciet_heater_v_2_0_steady_state.csv")
            .unwrap();

        wtr.write_record(&["time_seconds",
            "heater_power_kilowatts",
            "therminol_temperature_celsius",
            "shell_temperature_celsius",
            "timestep_seconds",])
            .unwrap();

        let mut time_wtr = Writer::from_path(
            "matrix_solve_trial_ciet_heater_v_2_0_steady_state_time.csv")
            .unwrap();

        time_wtr.write_record(&["loop_calculation_time_nanoseconds",
            "mutex_lock_time_ns",
            "node_connection_time_ns",
            "data_record_time_ns",
            "timestep_advance_time_ns",])
            .unwrap();

        let mut temp_profile_wtr = Writer::from_path(
            "matrix_solve_trial_ciet_heater_v_2_0_steady_state_temp_profile.csv")
            .unwrap();

        // this is code for writing the array of required temperatures
        {

            // I want the mid node length of this temperature

            let node_length: Length = heated_length/number_of_nodes as f64;

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


        // heater 2.0 while loop
        while current_time_simulation_time <= *max_time_ptr.deref(){

            let arc_mutex_lock_start = SystemTime::now();
            // obtain arc mutex locks 

            let mut back_cv_ptr_in_loop = 
            back_cv_ptr_clone_for_loop.lock().unwrap();
            let mut front_cv_ptr_in_loop = 
            front_cv_ptr_clone_for_loop.lock().unwrap();
            let mut fluid_temp_vec_ptr_in_loop
            = fluid_temp_vec_ptr_for_loop.lock().unwrap();
            let mut ambient_air_temp_at_present_timestep_ptr_in_loop
            = ambient_air_temp_at_present_timestep_ptr_for_loop.lock().unwrap();
            let mut q_fraction_ptr_in_loop 
            = q_fraction_ptr_for_loop.lock().unwrap();
            let mut fluid_vol_fraction_ptr_in_loop 
            = fluid_vol_fraction_ptr_for_loop.lock().unwrap();
            let mut solid_vol_fraction_ptr_in_loop 
            = solid_vol_fraction_ptr_for_loop.lock().unwrap();
            let mut steel_temperature_array_present_timestep_ptr_in_loop 
            = steel_temperature_array_present_timestep_ptr_for_loop.lock().unwrap();

            
            // time for arc mutex derefs
            let arc_mutex_lock_elapsed_ns = 
            arc_mutex_lock_start.elapsed().unwrap().as_nanos();
            
            // next, obtain average temperature 
            // average by volume for both fluid and solid vec
            //

            let node_connection_start = SystemTime::now();
            let mut fluid_temp_kelvin_times_vol_frac_average: 
            Vec<f64> = vec![];
            
            for (index,fluid_temp) in fluid_temp_vec_ptr_in_loop.iter().enumerate() {

                let fluid_temp_times_vol_frac = 
                fluid_vol_fraction_ptr_in_loop[index] *
                (*fluid_temp);

                fluid_temp_kelvin_times_vol_frac_average.push(
                    fluid_temp_times_vol_frac.get::<kelvin>());

            }

            // find average in kelvin, convert back to correct units 

            let fluid_avg_temp_kelvin: f64 = 
            fluid_temp_kelvin_times_vol_frac_average.iter().sum();

            let fluid_avg_temp = ThermodynamicTemperature::new::
                <kelvin>(fluid_avg_temp_kelvin);

            // we do the same for solid temp
            //
            //

            let mut steel_temp_kelvin_times_vol_frac_average: 
            Vec<f64> = vec![];
            
            for (index,steel_temp) in 
                steel_temperature_array_present_timestep_ptr_in_loop.iter().enumerate() {

                let steel_temp_times_vol_frac = 
                solid_vol_fraction_ptr_in_loop[index] *
                (*steel_temp);

                steel_temp_kelvin_times_vol_frac_average.push(
                    steel_temp_times_vol_frac.get::<kelvin>());

            }
            // heater 2.0
            // find average in kelvin, convert back to correct units 

            let steel_avg_temp_kelvin: f64 = 
            steel_temp_kelvin_times_vol_frac_average.iter().sum();

            let steel_avg_temp = ThermodynamicTemperature::new::
                <kelvin>(steel_avg_temp_kelvin);

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
                fluid_avg_temp,
                atmospheric_pressure);

            let therminol_steel_conductance_interaction: HeatTransferInteractionType
            = HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(
                    (steel, 
                    radial_thickness_thermal_conduction,
                    steel_avg_temp,
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
            get_thermal_conductance_based_on_interaction(
                fluid_avg_temp,
                steel_avg_temp,
                atmospheric_pressure,
                atmospheric_pressure,
                therminol_steel_conductance_interaction,
            ).unwrap();

            // now, create conductance vector  for heater v2.0

            let mut steel_therminol_conductance_vector: Array1<ThermalConductance> = 
            Array::zeros(number_of_nodes);

            steel_therminol_conductance_vector.fill(therminol_steel_nodal_thermal_conductance);

            // now we need conductance for steel and air, the procedure 
            // is quite similar
            let radial_thickness_thermal_conduction = 0.5*(od - id);
            let radial_thickness_thermal_conduction: 
            RadialCylindricalThicknessThermalConduction = 
            radial_thickness_thermal_conduction.into();

            let outer_diameter_thermal_conduction: OuterDiameterThermalConduction
            = od.into();

            let h_to_air: HeatTransfer = 
            HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);


            let steel_air_conductance_interaction: HeatTransferInteractionType
            = HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidOutside(
                    (steel, 
                    radial_thickness_thermal_conduction,
                    steel_avg_temp,
                    atmospheric_pressure),
                    (h_to_air,
                    outer_diameter_thermal_conduction,
                    node_length.clone().into())
            );

            let steel_air_nodal_thermal_conductance: ThermalConductance = 
            get_thermal_conductance_based_on_interaction(
                fluid_avg_temp,
                steel_avg_temp,
                atmospheric_pressure,
                atmospheric_pressure,
                steel_air_conductance_interaction,
            ).unwrap();

            // now, create conductance vector 

            let mut steel_air_conductance_vector: Array1<ThermalConductance> = 
            Array::zeros(number_of_nodes);

            steel_air_conductance_vector.fill(steel_air_nodal_thermal_conductance);


            // next is rho_cp vector for the fluid 
            // it has to be calculated based on each node temperature 
            // and should not be like the conductance bit
            //
            // hopefully it isn't too computationally expensive
            //

            let mut fluid_rho_cp_array: Array1<VolumetricHeatCapacity> 
            = Array::zeros(number_of_nodes);

            for (index,fluid_temp) in fluid_temp_vec_ptr_in_loop.iter().enumerate() {

                let fluid_nodal_rho_cp: VolumetricHeatCapacity = try_get_rho_cp(
                    therminol,
                    *fluid_temp,
                    atmospheric_pressure
                ).unwrap();

                fluid_rho_cp_array[index] = fluid_nodal_rho_cp;

            }

            // same for solid 
            let mut solid_rho_cp_array: Array1<VolumetricHeatCapacity> 
            = Array::zeros(number_of_nodes);

            for (index,solid_temp) in 
                steel_temperature_array_present_timestep_ptr_in_loop.iter().enumerate() {

                let solid_nodal_rho_cp: VolumetricHeatCapacity = try_get_rho_cp(
                    steel,
                    *solid_temp,
                    atmospheric_pressure
                ).unwrap();

                solid_rho_cp_array[index] = solid_nodal_rho_cp;

            }


            // now I'm going to manually make enthalpy inflows and 
            // outflows 

            let enthalpy_inflow_in_back_cv: Power 
            = therminol_mass_flowrate * try_get_h(
                therminol,
                inlet_temperature,
                atmospheric_pressure,
            ).unwrap();

            let enthalpy_outflow_in_front_cv: Power 
            = therminol_mass_flowrate * 
            front_cv_ptr_in_loop.deref().
            current_timestep_control_volume_specific_enthalpy;

            // add these to the power vectors inside each cv 

            back_cv_ptr_in_loop.deref_mut().rate_enthalpy_change_vector
            .push(enthalpy_inflow_in_back_cv);

            front_cv_ptr_in_loop.deref_mut().rate_enthalpy_change_vector
            .push(-enthalpy_outflow_in_front_cv);

            let node_connection_end_ns = 
            node_connection_start.elapsed().unwrap().as_nanos();

            // auto calculate timestep (not done)
            //
            // and also writing temperature profile 
            let data_recording_time_start = SystemTime::now();


            // record current fluid temperature profile
            {
                let mut temp_profile_data_vec: Vec<String> = vec![];

                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                // next simulation time string 
                let elapsed_calc_time_seconds_string = 
                calculation_time_elapsed.elapsed().unwrap().as_secs().to_string();
                temp_profile_data_vec.push(current_time_string);
                temp_profile_data_vec.push(elapsed_calc_time_seconds_string);

                for node_temp in steel_temperature_array_present_timestep_ptr_in_loop.deref_mut().iter() {

                    let node_temp_deg_c: f64 = 
                    node_temp.get::<degree_celsius>();

                    let node_temp_c_string: String = 
                    node_temp_deg_c.to_string();

                    temp_profile_data_vec.push(node_temp_c_string);
                }

                temp_profile_wtr.write_record(&temp_profile_data_vec).unwrap();
            }

            // outlet temperature profile 
            {
                // csv data writing

                let therminol_outlet_temp:ThermodynamicTemperature = 
                fluid_temp_vec_ptr_in_loop.deref()[number_of_nodes-1];

                let therminol_outlet_temp_string = 
                therminol_outlet_temp.get::<degree_celsius>().to_string();

                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                let heater_power_kilowatt_string = 
                heater_steady_state_power.get::<kilowatt>().to_string();

                // for st 11 

                // code block for recording inlet, shell and outlet 
                // temperatures
                //
                // now for shell temperatures, we are going to assume that 
                // ST-11 is used. 
                //
                // ST-11 is the thermocouple measuring surface temperature 
                // roughly 19 inches from the bottom of the heater 
                // The entire heated length excluding heater top and 
                // bottom heads is about 64 inches 
                //
                // So 19/64 is about 0.30 of the way through
                let st_11_length: Length = Length::new::<inch>(19_f64);


                // now I want to find out which node it is,
                // so i need the node length first 
                //
                
                let node_length: Length = heated_length/number_of_nodes as f64;

                // then use st_11 divide by node length 

                let st_11_rough_node_number: Ratio = st_11_length / node_length;

                // now, st_11 is about 19 inches, out of 64, and we have 
                // 8 equal nodes, each node is 
                // 12.5% of the heated length
                //
                // so this is about 30% of the way through
                //
                // so this is node three. 
                //
                // if we take st_11_length/node_length 
                // we would get about 2.375 for this ratio 
                //
                // we need to round up to get 3 
                // but the third node is the 2nd index in the matrix 
                // because the index starts from zero
                //
                //
                // so round it up and then minus 1 
                // most of the time, round down is ok, but rounding up 
                // makes more logical sense given this derivation

                let st_11_node_number: usize = 
                st_11_rough_node_number.get::<ratio>().ceil() as usize;

                let st_11_index_number: usize = st_11_node_number - 1;

                // now that we got the index number, we can get the 
                // outer surface temperature 

                let st_11_node_temp: ThermodynamicTemperature = 
                steel_temperature_array_present_timestep_ptr_in_loop
                .deref()[st_11_index_number];

                let shell_celsius_string = 
                st_11_node_temp.get::<degree_celsius>().to_string();
                
                // timestep in seconds
                //

                let timestep_string = 
                timestep.get::<second>().to_string();


                wtr.write_record(&[current_time_string,
                    heater_power_kilowatt_string,
                    therminol_outlet_temp_string,
                    shell_celsius_string,
                    timestep_string])
                    .unwrap();

            }


            let data_recording_time_end_ns = 
            data_recording_time_start.elapsed().unwrap().as_nanos();
            // now we are ready to start

            let mut fluid_temp_vec_last_timestep = 
            fluid_temp_vec_ptr_in_loop.deref().clone();

            let timestep_advance_start = 
            SystemTime::now();
            let _new_therminol_temperature_vec = 
            _advance_timestep_fluid_node_array_pipe_high_peclet_number(
                back_cv_ptr_in_loop.deref_mut(),
                front_cv_ptr_in_loop.deref_mut(),
                number_of_nodes,
                timestep,
                total_fluid_volume,
                Power::new::<watt>(0.0),
                steel_temperature_array_present_timestep_ptr_in_loop.deref_mut(),
                &mut steel_therminol_conductance_vector,
                fluid_temp_vec_ptr_in_loop.deref_mut(),
                therminol_mass_flowrate,
                fluid_vol_fraction_ptr_in_loop.deref_mut(),
                &mut fluid_rho_cp_array,
                q_fraction_ptr_in_loop.deref_mut(),
            ).unwrap();

            let _new_steel_temperature_vec = 
            shell_solid_node::_advance_timestep_solid_cylindrical_shell_node_no_axial_conduction(
                number_of_nodes,
                timestep,
                total_steel_volume,
                heater_steady_state_power,
                &mut fluid_temp_vec_last_timestep,
                &mut steel_therminol_conductance_vector,
                ambient_air_temp_at_present_timestep_ptr_in_loop.deref_mut(),
                &mut steel_air_conductance_vector,
                steel_temperature_array_present_timestep_ptr_in_loop.deref_mut(),
                solid_vol_fraction_ptr_in_loop.deref_mut(),
                &mut solid_rho_cp_array,
                &mut q_fraction_ptr_in_loop
            ).unwrap();

            
            let timestep_advance_end_ns = 
            timestep_advance_start.elapsed().unwrap().as_nanos();

            // write timestep diagnostics

            let total_time_ns = 
            arc_mutex_lock_elapsed_ns 
            + node_connection_end_ns 
            + data_recording_time_end_ns 
            + timestep_advance_end_ns;

            time_wtr.write_record(&[total_time_ns.to_string(),
                arc_mutex_lock_elapsed_ns.to_string(),
                node_connection_end_ns.to_string(),
                data_recording_time_end_ns.to_string(),
                timestep_advance_end_ns.to_string()])
                .unwrap();
            // timestep addition to current simulation time
            current_time_simulation_time += timestep;
        }
        // with csvs being written,
        // use cargo watch -x test --ignore '*.csv'
        time_wtr.flush().unwrap();
        temp_profile_wtr.flush().unwrap();
        wtr.flush().unwrap();
    };

    let main_thread_handle = thread::spawn(calculation_loop);

    main_thread_handle.join().unwrap();


    return ();

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
    let mu: DynamicViscosity = try_get_mu_viscosity(therminol,
        therminol_temperature,
        pressure).unwrap();

    let k: ThermalConductivity = try_get_kappa_thermal_conductivity(
        therminol,
        therminol_temperature,
        pressure).unwrap();


    let reynolds: Ratio = _ciet_heater_v_2_0_reynolds_nunber(
        mass_flowrate, mu);

    let prandtl: Ratio = try_get_prandtl(
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
