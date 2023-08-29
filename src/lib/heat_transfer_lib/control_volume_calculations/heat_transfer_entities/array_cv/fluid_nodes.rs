use std::{sync::{Arc, Mutex}, time::SystemTime, ops::Deref, thread};

use ndarray::*;
use ndarray_linalg::{*, error::LinalgError};
use uom::{si::{f64::*, power::{watt, kilowatt}, length::{meter, centimeter}, thermodynamic_temperature::{degree_celsius, kelvin}, pressure::atmosphere, area::square_meter, mass_rate::kilogram_per_second, time::second, heat_transfer::watt_per_square_meter_kelvin}, num_traits::Zero};

use crate::heat_transfer_lib::{thermophysical_properties::{specific_enthalpy::specific_enthalpy, LiquidMaterial, SolidMaterial}};
use crate::heat_transfer_lib::thermophysical_properties::Material;
use crate::heat_transfer_lib::control_volume_calculations::{heat_transfer_interactions::HeatTransferInteractionType};
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::{SingleCVNode, array_cv::calculation::solve_conductance_matrix_power_vector, HeatTransferEntity, RadialCylindricalThicknessThermalConduction, InnerDiameterThermalConduction};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::CVType::*;

/// for high peclet number flows, we can advance timestep without 
/// considering axial conduction
///
/// However, this fluid node array will be connected to a shell of some 
/// sort, like a metallic pipe 
///
/// There will be a conductance array which contains information of how 
/// much thermal conductance connects the inner nodes to the outer 
/// metallic pipe array
///
/// [solid]
/// T_back            T[0]           T[1]          T[n-1]         T_front 
///
/// ----------------fluid solid boundary with thermal resistance ---------
///
/// [fluid]
/// T_back            T[0]           T[1]          T[n-1]         T_front 
///
/// At the back and front node, there will be a single_cv which is able 
/// to link up to other cvs
///
/// You can also add heat volumetrically to the fluid as if it were 
/// generating heat, but you can also set it to zero
///
///
pub (in crate) 
fn advance_timestep_fluid_node_array_pipe_high_peclet_number(
    back_single_cv: &mut SingleCVNode,
    front_single_cv: &mut SingleCVNode,
    number_of_nodes: usize,
    dt: Time,
    total_volume: Volume,
    q: Power,
    last_timestep_temperature_solid: &mut Array1<ThermodynamicTemperature>,
    solid_fluid_conductance_array: &mut Array1<ThermalConductance>,
    last_timestep_temperature_fluid: &mut Array1<ThermodynamicTemperature>,
    mass_flowrate: MassRate,
    volume_fraction_array: &mut Array1<f64>,
    rho_cp: &mut Array1<VolumetricHeatCapacity>,
    q_fraction: &mut Array1<f64>)
-> Result<Array1<ThermodynamicTemperature>,error::LinalgError>{

    // First things first, we need to set up 
    // how the CV interacts with the internal array
    // here is heat added to CV

    let back_cv_rate_enthalpy_change_vector: Vec<Power> = 
    back_single_cv.rate_enthalpy_change_vector.clone();

    let front_cv_rate_enthalpy_change_vector: Vec<Power> = 
    front_single_cv.rate_enthalpy_change_vector.clone();

    // compute power source for back node

    let mut total_enthalpy_rate_change_back_node = 
    Power::new::<watt>(0.0);

    for enthalpy_chg_rate in 
        back_cv_rate_enthalpy_change_vector.iter() {

            total_enthalpy_rate_change_back_node += *enthalpy_chg_rate;
        }

    // then the front node,

    let mut total_enthalpy_rate_change_front_node = 
    Power::new::<watt>(0.0);

    for enthalpy_chg_rate in 
        front_cv_rate_enthalpy_change_vector.iter() {

            total_enthalpy_rate_change_front_node += *enthalpy_chg_rate;
        }

    // this front and back nodes will be an extra term added to the 
    // heat source vector S
    //
    // The old fluid temperature will need to be used to calculate 
    // new specific enthalpy for the system

    let mut temperature_vector: Array1<ThermodynamicTemperature> = 
    last_timestep_temperature_fluid.map(
        |temp_kelvin_ptr: &ThermodynamicTemperature| {

            return *temp_kelvin_ptr;
        }

    );

    // now let's start calculation 
    //
    // there will always be at least 2 nodes

    if number_of_nodes <= 1 {
        return Err(LinalgError::Shape(
            ShapeError::from_kind(
                ErrorKind::OutOfBounds
            )));
    } else if number_of_nodes > 1 {

        let mut coefficient_matrix: Array2<ThermalConductance> = 
        Array::zeros((number_of_nodes, number_of_nodes));

        let mut power_source_vector: 
        Array1<Power> = Array::zeros(number_of_nodes);

        // ascertain if we have forward flow 

        let forward_flow: bool = mass_flowrate.ge(&MassRate::zero());

        // back node calculation (first node)
        {
            // for the first node, also called the back node
            // energy balance is: 
            // m c_p dT/dt = -H (T - T_solid) - m_flow h_fluid(T) 
            // + m_flow h_fluid(adjacent T) + q
            // of all these terms, only the m cp dT/dt term and HT 
            // 
            // We separate this out to get:
            //
            // m cp T / dt + HT = 
            //
            // HT_solid - m_flow h_fluid(T_old) 
            // + m_flow h_fluid(adjacent T_old) + m cp / dt (Told)
            // + q
            //
            // belong in the M matrix, the rest belong in S
            coefficient_matrix[[0,0]] = volume_fraction_array[0] * rho_cp[0] 
            * total_volume / dt + solid_fluid_conductance_array[0];

            // the first part of the source term deals with 
            // the flow direction independent terms

            let h_fluid_last_timestep: AvailableEnergy = 
            back_single_cv.current_timestep_control_volume_specific_enthalpy;

            // now this makes the scheme semi implicit, and we should then 
            // treat the scheme as explicit

            power_source_vector[0] = solid_fluid_conductance_array[0] *
                last_timestep_temperature_solid[0] 
                - mass_flowrate * h_fluid_last_timestep 
                + last_timestep_temperature_fluid[0] * total_volume * 
                volume_fraction_array[0] * rho_cp[0] / dt 
                + q * q_fraction[0] 
                + total_enthalpy_rate_change_back_node ;

            // the next part deals with the inflow
            // m_flow h_fluid(adjacent T_old)
            //
            // now if the advection interaction is done correctly, 
            //
            // (advection) ----- (back cv) --------> fwd
            //
            // then in a frontal flow condition, the enthalpy flows in 
            // would already have been accounted for
            //
            // but in the case of backflow, then fluid from the node
            // in front will flow into this fluid node 
            // that is node 1 

            // so if mass flowrate is <= 0 , then we will calculate 
            // backflow conditions

            if !forward_flow {
                // first, get enthalpy of the node in front 

                let enthalpy_of_adjacent_node_to_the_front: AvailableEnergy = 
                specific_enthalpy(
                    back_single_cv.material_control_volume,
                    last_timestep_temperature_fluid[1],
                    back_single_cv.pressure_control_volume).unwrap();

                // now if mass flowrate is less than zero, then 
                // we receive enthalpy from the front cv 
                //
                // But we need to subtract the negative mass flow if 
                // that makes sense, or at least make it absolute
                // - (-m) * h = m * h
                //

                power_source_vector[0] += 
                mass_flowrate.abs() * enthalpy_of_adjacent_node_to_the_front;



            }
        }

        // bulk node calculations 
        if number_of_nodes > 2 {
            // loop over all nodes from 1 to n-2 (n-1 is not included)
            for i in 1..number_of_nodes-1 {
                // the coefficient matrix is pretty much the same, 
                // we only consider solid fluid conduction, and the 
                // thermal inertia terms

                coefficient_matrix[[i,i]] = volume_fraction_array[i] * rho_cp[i] 
                    * total_volume / dt + solid_fluid_conductance_array[i];

                // we also consider outflow using previous timestep 
                // temperature, 
                // assume back cv and front cv material are the same

                let h_fluid_last_timestep: AvailableEnergy = 
                specific_enthalpy(
                    back_single_cv.material_control_volume,
                    last_timestep_temperature_fluid[i],
                    back_single_cv.pressure_control_volume).unwrap();

                // basically, all the power terms remain 
                power_source_vector[i] = solid_fluid_conductance_array[i] *
                    last_timestep_temperature_solid[i] 
                    - mass_flowrate * h_fluid_last_timestep 
                    + last_timestep_temperature_fluid[i] * total_volume * 
                    volume_fraction_array[i] * rho_cp[i] / dt 
                    + q * q_fraction[i];

                // account for enthalpy inflow

                if forward_flow {

                    // enthalpy must be based on the the cv at i-1

                    let h_fluid_adjacent_node: AvailableEnergy = 
                    specific_enthalpy(
                        back_single_cv.material_control_volume,
                        last_timestep_temperature_fluid[i-1],
                        back_single_cv.pressure_control_volume).unwrap();

                    
                    power_source_vector[i] += 
                    h_fluid_adjacent_node * mass_flowrate.abs();

                } else {

                    // enthalpy must be based on cv at i+1
                    let h_fluid_adjacent_node: AvailableEnergy = 
                    specific_enthalpy(
                        back_single_cv.material_control_volume,
                        last_timestep_temperature_fluid[i+1],
                        back_single_cv.pressure_control_volume).unwrap();

                    
                    power_source_vector[i] += 
                    h_fluid_adjacent_node * mass_flowrate.abs();

                }

            }
        }

        // front node (last node) calculation 
        {
            // we still follow the same guideline
            // m cp T / dt + HT = 
            //
            // HT_solid - m_flow h_fluid(T_old) 
            // + m_flow h_fluid(adjacent T_old) + m cp / dt (Told)
            // + q
            let i = number_of_nodes-1;

            coefficient_matrix[[i,i]] = volume_fraction_array[i] * rho_cp[i] 
            * total_volume / dt + solid_fluid_conductance_array[i];
            // the first part of the source term deals with 
            // the flow direction independent terms

            let h_fluid_last_timestep: AvailableEnergy = 
            front_single_cv.current_timestep_control_volume_specific_enthalpy;
            // now this makes the scheme semi implicit, and we should then 
            // treat the scheme as explicit

            power_source_vector[i] = solid_fluid_conductance_array[i] *
                last_timestep_temperature_solid[i] 
                - mass_flowrate * h_fluid_last_timestep 
                + last_timestep_temperature_fluid[i] * total_volume * 
                volume_fraction_array[i] * rho_cp[i] / dt 
                + q * q_fraction[i] 
                + total_enthalpy_rate_change_front_node ;

            // now advection causing heat transfer between the frontal 
            // node and some other single cv has already been accounted 
            // for in total_enthalpy_rate_change_front_node 
            // If flow is forward facing though, then enthalpy will 
            // be coming in from the i-1 th node

            if forward_flow {
                // first, get enthalpy of the adjacent node to the 
                // back

                let enthalpy_of_adjacent_node_to_the_rear: AvailableEnergy = 
                specific_enthalpy(
                    back_single_cv.material_control_volume,
                    last_timestep_temperature_fluid[i-1],
                    back_single_cv.pressure_control_volume).unwrap();

                // now if mass flowrate is less than zero, then 
                // we receive enthalpy from the front cv 
                //
                // But we need to subtract the negative mass flow if 
                // that makes sense, or at least make it absolute
                // - (-m) * h = m * h
                //

                power_source_vector[i] += 
                mass_flowrate.abs() * enthalpy_of_adjacent_node_to_the_rear;

                // if there's backflow, there's no further enthalpy 
                // inflow 

            }
        }

        // solve for new temperature 

        temperature_vector = 
            solve_conductance_matrix_power_vector(
                coefficient_matrix,power_source_vector)?;
    } 
    // update the single cvs at the front and back with new enthalpies 

    // Todo: probably need to synchronise error types in future
    let back_node_enthalpy_next_timestep: AvailableEnergy = 
    specific_enthalpy(
        back_single_cv.material_control_volume,
        temperature_vector[0],
        back_single_cv.pressure_control_volume).unwrap();

    back_single_cv.current_timestep_control_volume_specific_enthalpy 
    = back_node_enthalpy_next_timestep;

    let front_node_enthalpy_next_timestep: AvailableEnergy = 
    specific_enthalpy(
        front_single_cv.material_control_volume,
        temperature_vector[number_of_nodes-1],
        front_single_cv.pressure_control_volume).unwrap();

    front_single_cv.current_timestep_control_volume_specific_enthalpy 
    = front_node_enthalpy_next_timestep;

    // let's also update the previous temperature vector 

    *last_timestep_temperature_fluid = temperature_vector.mapv(
        |temperature_value| {
            return temperature_value;
        }
    );

    // set liquid cv mass 
    // probably also need to update error types in future
    back_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
    back_single_cv.rate_enthalpy_change_vector.clear();
    back_single_cv.max_timestep_vector.clear();

    front_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
    front_single_cv.rate_enthalpy_change_vector.clear();
    front_single_cv.max_timestep_vector.clear();

    return Ok(temperature_vector);

}

#[test]
fn fluid_node_calculation_initial_test(){


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
    let timestep_placeholder = Time::new::<second>(0.01);
    let total_volume = flow_area * heated_length;


    let mut initial_temperature_array: Array1<ThermodynamicTemperature> = 
    Array::default(number_of_nodes);
    initial_temperature_array.fill(initial_temperature);

    let fluid_temperature_array = initial_temperature_array.clone();

    let mut steel_temperature_array: Array1<ThermodynamicTemperature> = 
    Array::default(number_of_nodes);

    steel_temperature_array.fill(ambient_air_temp);




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
    let steel_temp_at_present_timestep_ptr = Arc::new(Mutex::new(
        steel_temperature_array
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

    // clone the ptrs for the loop
    let back_cv_ptr_clone_for_loop = back_cv_ptr.clone();
    let front_cv_ptr_clone_for_loop = front_cv_ptr.clone();
    let fluid_temp_vec_ptr_for_loop = fluid_temp_at_present_timestep_ptr.clone();
    let steel_temp_at_present_timestep_ptr_for_loop = 
    steel_temp_at_present_timestep_ptr.clone();
    let q_fraction_ptr_for_loop = q_fraction_ptr.clone();
    let fluid_vol_fraction_ptr_for_loop = fluid_vol_fraction_ptr.clone();
    let solid_vol_fraction_ptr_for_loop = solid_vol_fraction_ptr.clone();

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

            let back_cv_ptr_in_loop = 
            back_cv_ptr_clone_for_loop.lock().unwrap();
            let front_cv_ptr_in_loop = 
            front_cv_ptr_clone_for_loop.lock().unwrap();
            let fluid_temp_vec_ptr_in_loop
            = fluid_temp_vec_ptr_for_loop.lock().unwrap();
            let steel_temp_at_present_timestep_ptr_in_loop
            = steel_temp_at_present_timestep_ptr_for_loop.lock().unwrap();
            let q_fraction_ptr_in_loop 
            = q_fraction_ptr_for_loop.lock().unwrap();
            let fluid_vol_fraction_ptr_in_loop 
            = fluid_vol_fraction_ptr_for_loop.lock().unwrap();
            let solid_vol_fraction_ptr_in_loop 
            = solid_vol_fraction_ptr_for_loop.lock().unwrap();

            // advance timestep
            current_time_simulation_time += timestep;
            
            // time for arc mutex derefs
            let arc_mutex_lock_end = SystemTime::now();

            let arc_mutex_lock_elapsed_ms = 
            arc_mutex_lock_start.elapsed().unwrap().as_millis()
            - arc_mutex_lock_end.elapsed().unwrap().as_millis();
            
            // next, obtain average temperature 
            // average by volume for both fluid and solid vec
            //

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
                steel_temp_at_present_timestep_ptr_in_loop.iter().enumerate() {

                let steel_temp_times_vol_frac = 
                solid_vol_fraction_ptr_in_loop[index] *
                (*steel_temp);

                steel_temp_kelvin_times_vol_frac_average.push(
                    steel_temp_times_vol_frac.get::<kelvin>());

            }

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

            let h: HeatTransfer = 
            HeatTransfer::new::<watt_per_square_meter_kelvin>(35.0);

            let conductance_interaction: HeatTransferInteractionType
            = HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(
                    (steel, 
                    radial_thickness_thermal_conduction,
                    steel_avg_temp,
                    atmospheric_pressure),
                    (h,
                    inner_diameter_thermal_conduction,
                    node_length.clone().into())
            );

            // now based on conductance interaction, 
            // we can obtain thermal conductance, the temperatures 
            // and pressures don't really matter


            
        }
    };

    let main_thread_handle = thread::spawn(calculation_loop);

    main_thread_handle.join().unwrap();


    return ();

}
