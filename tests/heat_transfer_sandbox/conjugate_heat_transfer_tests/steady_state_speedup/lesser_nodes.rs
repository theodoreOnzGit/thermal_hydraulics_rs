use std::ops::{DerefMut, Deref};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::SystemTime;

use csv::Writer;
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
    SurfaceArea, SingleCVNode, CVType, BCType, InnerDiameterThermalConduction, OuterDiameterThermalConduction, RadialCylindricalThicknessThermalConduction, CylinderLengthThermalConduction};
use thermal_hydraulics_rs::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::CVType::SingleCV;
use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::density::density;



use uom::si::angle::radian;
use uom::si::angular_velocity::radian_per_second;
use uom::si::area::{square_centimeter, square_meter};
use uom::si::f64::*;
use uom::si::length::{centimeter, meter, inch};
use uom::si::mass_rate::kilogram_per_second;
use uom::si::power::{watt, kilowatt};
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

use rayon::prelude::*;

/// In this test, we have a nodalised representation of the 
/// this is 8 nodes in the axial direction and 2 nodes for metal 
/// in the radial direction
///
/// This is for heater v2.0
///
/// I realised that we needed a major speedup because in debug mode 
/// a 2.3ms timestep was taking ~ 602 ms to calculate
///
/// This was unacceptable
///
/// Now, I introduced partial parallelism in this test, because 
/// it was quite cumbersome to do. This is of course experimental 
///
/// I spawned 3 threads to connect the fluid nodes to the inner nodes 
/// and the time per timestep in debug mode went from 602 ms to about 
/// 562 ms
///
/// When 4 threads were used, it reduced to about 560 ms 
/// this was on an intel i7-10875H with 32GB ram, linux system
/// In this case, the thread::spawn is expected to use some overhead 
/// and therefore, it will be somewhat slow
///
/// we can either minimise the overhead by spawning a rayon threadpool
/// (ie spawn a persistent thread ONCE at process start, not in the 
/// while loop)
/// or just live with it. Which means minimise thread spawning in Rust,
/// maybe three or four threads, at most in each while loop
///
/// The former means we don't need the extra dependency
/// so we can just rely on the std library. It's also easier to understand,
/// but much more cumbersome to code 
///
/// However, with thread spawning, one must be aware of how to 
/// handle this 
///
/// This is because mutable vectors can only mutate one element at a time.
/// The best is to use Arc Mutex pointers and clone each individual heat 
/// transfer entity, and make at least two of those pointers in one go. 
/// One pointer gets absorbed in the thread, and is destroyed because it 
/// is moved there, the other remains outside the loop and should be 
/// derferenced at the last part to update the vectors
///
/// Now for this, I only parallelised the fluid node to first layer 
/// of shell connection 
/// 
/// (ambient air)
///
/// | 
/// | 
/// | 
/// 
/// (steel outer shell)  ---- (heater power)
///
/// | 
/// | 
/// | 
///
/// (steel inner shell) ---- (heater power)
///
/// | 
/// |  <---- parallelised, rest is serial calculation
/// | 
///
/// (fluid node)  --- (subsequent nodes, advection)
///
///
///
/// so one out of six interactions is parallelised, it saved 40ms out 
/// of roughly 400 ms of calculation time, roughly 10%.
///
/// If we extrapolate all this, we can expect to save perhaps 50% of the 
/// simulation time here.
///
/// However, for individual heat transfer elements, it's kind of cumbersome 
/// to do things like this. However, you still got to connect roughly 
/// 30 or so components together in CIET, and using thread spawn for 
/// parallel joining will be useful in that area.
///
/// If we parallelise all the node connections with thread spawn, maybe 
/// we can expect 50% time saving? So it's good but not a cure all.
/// Nevertheless, tool is indispensable, an we are still using slow cloning 
/// processes, and not using Rayon. 
/// 
///
#[test]
//#[ignore = "data collected"]
pub fn ciet_heater_v_2_0_test_steady_state_v_1_2_three_nodes(){


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
    let number_of_nodes: usize = 3;
    let ambient_air_temp = ThermodynamicTemperature::new::<
        degree_celsius>(21.67);

    let heater_steady_state_power = Power::new::<kilowatt>(8.0);

    
    

    let therminol_mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);

    // construct the objects,
    // I'm going to use a function 

    fn construct_heated_section_fluid_nodes(therminol: Material,
        cross_sectional_area: Area,
        heated_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let mut fluid_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = heated_length/number_of_nodes as f64;

        for _index in 0..number_of_nodes {
            let therminol_node: HeatTransferEntity = 
            SingleCVNode::new_odd_shaped_pipe(
                node_length,
                cross_sectional_area,
                therminol,
                initial_temperature,
                pressure,
            ).unwrap();

            fluid_node_vec.push(therminol_node);
        }

        return fluid_node_vec;
    }

    let fluid_node_vec: Vec<HeatTransferEntity> 
    = construct_heated_section_fluid_nodes(
        therminol,
        flow_area,
        heated_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);
    
    // then construct two layers of steel shells
    
    fn construct_steel_shell_nodes(steel: Material,
        id: Length, 
        od: Length,
        heated_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let id: InnerDiameterThermalConduction = id.into();
        let od: OuterDiameterThermalConduction = od.into();

        let mut steel_shell_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = heated_length/number_of_nodes as f64;

        for _index in 0..number_of_nodes {

            let steel_shell_node = SingleCVNode::new_cylindrical_shell(
                node_length,
                id, od,
                steel, 
                initial_temperature,
                pressure,
            ).unwrap();

            steel_shell_node_vec.push(steel_shell_node);
        }

        return steel_shell_node_vec;
    }

    // inner layer of steel shell I will just assume 0.0392 m is the 
    // midway point
    // the inner node should be thicker anyway 

    let midway_point_steel_shell: Length = 
    Length::new::<meter>(0.0392);

    let steel_shell_inner_node_vec: Vec<HeatTransferEntity> = 
    construct_steel_shell_nodes(
        steel,
        id, midway_point_steel_shell, heated_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);

    let steel_shell_outer_node_vec: Vec<HeatTransferEntity> = 
    construct_steel_shell_nodes(
        steel,
        midway_point_steel_shell, od, heated_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);

    // now, let me make mutex locks and Arc pointers

    let fluid_node_vec_ptr = Arc::new(Mutex::new(
        fluid_node_vec
    ));

    let steel_shell_inner_node_vec_ptr = Arc::new(Mutex::new(
        steel_shell_inner_node_vec
    ));

    let steel_shell_outer_node_vec_ptr = Arc::new(Mutex::new(
        steel_shell_outer_node_vec
    ));
    






    // need two boundary conditions 

    let inlet_const_temp = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ThermodynamicTemperature::new::<degree_celsius>(79.12)
        ));

    let outlet_zero_heat_flux = HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatAddition(Power::new::<watt>(0.0))
    );

    let ambient_temperature_bc = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ambient_air_temp
        ));

    let inlet_const_temp_ptr = Arc::new(Mutex::new(
        inlet_const_temp
    ));

    let outlet_zero_heat_flux_ptr = Arc::new(Mutex::new(
        outlet_zero_heat_flux
    ));

    let ambient_air_temp_bc_ptr = Arc::new(Mutex::new(
        ambient_temperature_bc
    ));

    // the two types of HeatTransferInteractionType are 
    // advection and convection resistance
    //
    // 2007 square_centimeter
    // and 607 watt_per_square_meter_kelvin



    // timestep settings


    let max_time: Time = Time::new::<second>(200.0);
    let max_time_ptr = Arc::new(max_time);

    let calculation_time_elapsed = SystemTime::now();

    // this is the calculation loop
    let calculation_loop = move || {


        // csv writer, for post processing 

        // now for heater power, can swap between 

        let heater_power = heater_steady_state_power;


        let mut wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_three_axial_nodes.csv")
            .unwrap();

        wtr.write_record(&["time_seconds",
            "heater_power_kilowatts",
            "therminol_temperature_celsius",
            "shell_temperature_celsius",
            "auto_timestep_calculated_seconds",])
            .unwrap();

        let mut time_wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_time_three_axial_nodes.csv")
            .unwrap();

        time_wtr.write_record(&["loop_calculation_time_ms",
            "mutex_lock_frac",
            "node_connection_frac",
            "data_record_frac",
            "timestep_advance_frac",])
            .unwrap();

        // now i want a writer for temperature profile over n nodes 
        // it needs a simulation time, computation time elapsed, 
        // and a temperature for the outer surface temperature node 
        // for all nodes 

        let mut temp_profile_wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_temp_profile_three_axial_nodes.csv")
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


        let mut current_time_simulation_time = Time::new::<second>(0.0);

        let max_time_ptr_in_loop = max_time_ptr;
        // we are sampling at about 10 Hz
        // so the nyquist frequency is about 5 Hz 
        // this is because the highest frequency is about 3.66 Hz
        
        while current_time_simulation_time <= *max_time_ptr_in_loop {

            // timer for timekeeping purposes 

            let arc_mutex_lock_start = SystemTime::now();

            // calculation steps

            let mut inlet_const_temp_in_loop = 
            inlet_const_temp_ptr.lock().unwrap();
            let mut outlet_zero_heat_flux_in_loop = 
            outlet_zero_heat_flux_ptr.lock().unwrap();
            
            let mut ambient_air_temp_bc_in_loop = 
            ambient_air_temp_bc_ptr.lock().unwrap();

            let mut fluid_vec_in_loop = 
            fluid_node_vec_ptr.lock().unwrap();

            let mut steel_shell_inner_node_vec_in_loop = 
            steel_shell_inner_node_vec_ptr.lock().unwrap();

            let mut steel_shell_outer_node_vec_in_loop = 
            steel_shell_outer_node_vec_ptr.lock().unwrap();

            let arc_mutex_lock_end = SystemTime::now();

            let arc_mutex_lock_elapsed_ms = 
            arc_mutex_lock_start.elapsed().unwrap().as_millis()
            - arc_mutex_lock_end.elapsed().unwrap().as_millis();

            // we need to conenct a few things 
            //
            // to simplify, we also ignore axial conduction 
            // in all materials
            // 
            // in the radial direction:
            //
            // (ambient temp)
            // |
            // | --  convection/conduction
            // |
            // (outer steel shell nodes)
            // |
            // | -- conduction (cylindrical)
            // |
            // (inner steel shell nodes) 
            // | 
            // | -- conduction/convection
            // | 
            // (fluid nodes) 
            //
            // 
            // in the axial direction: 
            //
            // (inlet) --- (fluid nodes) --- (outlet)

            // calculate convection interactions
            // first, we need to connect the 
            //

            let node_connection_start = SystemTime::now();


            // radial heat transfer interactions
            // this one is fluid to inner steel cylinder
            for (index, fluid_node_heat_trf_entity_ptr) in 
                fluid_vec_in_loop.iter_mut().enumerate(){

                    // conduction material, properties

                    // radial thickness needs to be half of the 
                    // shell thickness, because it links to the shell 
                    // center
                    let radial_thickness: Length = 
                    (midway_point_steel_shell - id) *0.5;

                    let radial_thickness: RadialCylindricalThicknessThermalConduction
                    = radial_thickness.into();

                    let steel_inner_cylindrical_node_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature(
                        &mut steel_shell_inner_node_vec_in_loop[index]).unwrap();

                    let pressure = atmospheric_pressure;

                    // now need to get heat transfer coeff

                    let therminol_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature(
                        fluid_node_heat_trf_entity_ptr).unwrap();

                    let heat_trf_coeff: HeatTransfer = 
                    heat_transfer_coefficient_ciet_v_2_0(
                        therminol_mass_flowrate,
                        therminol_temp,
                        pressure,
                    );

                    let inner_diameter: InnerDiameterThermalConduction = 
                    id.clone().into();

                    let node_length: Length = 
                    heated_length/(number_of_nodes as f64);

                    let node_length: CylinderLengthThermalConduction = 
                    node_length.into();

                    // construct the interaction 

                    let interaction: HeatTransferInteractionType = 
                    HeatTransferInteractionType::CylindricalConductionConvectionLiquidInside
                        ((steel,radial_thickness,
                            steel_inner_cylindrical_node_temp,
                            pressure),
                            (heat_trf_coeff,
                                inner_diameter,
                                node_length));

                    // link the entities,
                    // this is the fluid to the inner shell
                    link_heat_transfer_entity(fluid_node_heat_trf_entity_ptr, 
                        &mut steel_shell_inner_node_vec_in_loop[index], 
                        interaction).unwrap();

                }

            // second, link inner shell to outer shell 
            //
            // both inner and outer shell will have power as well
            // roughly evenly distributed
            for (index, inner_shell_ptr) in steel_shell_inner_node_vec_in_loop.
                iter_mut().enumerate(){
                    let outer_shell_ptr: &mut HeatTransferEntity = 
                    &mut steel_shell_outer_node_vec_in_loop[index];


                    // remember, the inner diameter and outer diameter 
                    // of the shells are at the midpoints of both shells 
                    //
                    // so the radial thicknesses are halved.
                    let inner_radial_thickness: Length = 
                    0.5*(midway_point_steel_shell - id);

                    let outer_radial_thickness: Length = 
                    0.5*(od - midway_point_steel_shell);

                    let inner_radial_thickness: RadialCylindricalThicknessThermalConduction 
                    = inner_radial_thickness.into();

                    let outer_radial_thickness: RadialCylindricalThicknessThermalConduction 
                    = outer_radial_thickness.into();

                    // remember, the inner diameter and outer diameter 
                    // of the shells are at the midpoints of both shells 
                    //
                    // thats why you see all this subtraction

                    let id_mid_inner_shell: Length = 
                    id + inner_radial_thickness.into();

                    let id_mid_inner_shell: InnerDiameterThermalConduction 
                    = id_mid_inner_shell.into();

                    let od_mid_inner_shell: Length = 
                    od - outer_radial_thickness.into();

                    let od_mid_inner_shell: OuterDiameterThermalConduction 
                    = od_mid_inner_shell.into();

                    let node_length: Length = 
                    heated_length/(number_of_nodes as f64);

                    let node_length: CylinderLengthThermalConduction = 
                    node_length.into();

                    // create the interaction 

                    let interaction = HeatTransferInteractionType::
                        DualCylindricalThermalConductance(
                            (steel, inner_radial_thickness),
                            (steel, outer_radial_thickness),
                            (id_mid_inner_shell, od_mid_inner_shell, node_length),
                        );

                    // link them together
                    link_heat_transfer_entity(inner_shell_ptr, 
                        outer_shell_ptr, 
                        interaction).unwrap();

                    // each node would have roughly the same power 
                    // I should of course average it volumetrically 
                    // but I'm not going to
                    // 
                    // I'll divide it by the number of nodes, 
                    // for the number of nodes in each shell layer 
                    // and then divide by the number of layers
                    //
                    // I'll have to immutably borrow the fluid_vec_in_loop
                    // ptr because the borrowing rules make it hard 
                    // to borrow outer shell or inner shell
                    let node_heater_power: Power;
                    {
                        node_heater_power = 
                        heater_power / 2 as f64  
                        / fluid_vec_in_loop.len() as f64;
                    }

                    let mut electrical_heat_bc: HeatTransferEntity = 
                    BCType::new_const_heat_addition(node_heater_power);

                    let heat_addition_interaction = 
                    HeatTransferInteractionType::UserSpecifiedHeatAddition;

                    // link the power BC to the inner shell
                    link_heat_transfer_entity(inner_shell_ptr, 
                        &mut electrical_heat_bc, 
                        heat_addition_interaction).unwrap();

                    // link them together
                    link_heat_transfer_entity(outer_shell_ptr, 
                        &mut electrical_heat_bc, 
                        heat_addition_interaction).unwrap();

                }

            // third, link outer shell to ambient temperature

            for (index,outer_shell_ptr) in 
                steel_shell_outer_node_vec_in_loop.iter_mut().enumerate(){

                    // conduction material, properties

                    // radial thickness needs to be half of the 
                    // shell thickness, because it links to the shell 
                    // center
                    //
                    // (outer shell) ----- (outer shell surface) ---- (fluid)
                    // 
                    //          R_{shell} (half length)    R_{conv}
                    //
                    //
                    //  for convection, h = 20 W/(m^2 K)
                    //
                    let radial_thickness: Length = 
                    (od - midway_point_steel_shell) *0.5;

                    let radial_thickness: RadialCylindricalThicknessThermalConduction
                    = radial_thickness.into();

                    let steel_outer_cylindrical_node_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature(
                        &mut steel_shell_inner_node_vec_in_loop[index]).unwrap();

                    let pressure = atmospheric_pressure;

                    // now need to get heat transfer coeff
                    //
                    // 20 W/(m^2 K)

                    let heat_trf_coeff: HeatTransfer = 
                    HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);


                    let outer_diameter: OuterDiameterThermalConduction = 
                    od.clone().into();

                    let node_length: Length = 
                    heated_length/(number_of_nodes as f64);

                    let node_length: CylinderLengthThermalConduction = 
                    node_length.into();

                    // construct the interaction 

                    let interaction: HeatTransferInteractionType = 
                    HeatTransferInteractionType::CylindricalConductionConvectionLiquidOutside
                        ((steel,radial_thickness,
                            steel_outer_cylindrical_node_temp,
                            pressure),
                            (heat_trf_coeff,
                                outer_diameter,
                                node_length));

                    // link the entities,
                    // this is the fluid to the inner shell
                    link_heat_transfer_entity(&mut ambient_air_temp_bc_in_loop, 
                        outer_shell_ptr, 
                        interaction).unwrap();

                    


                }
            
            // fourth, link adjacent fluid nodes axially with advection 

            for index in 0..fluid_vec_in_loop.len()-1 {

                // (heater node i) ----- (heater node 2) ----- .... 
                //
                //
                //
                // to borrow two elements mutably is tricky, 
                // so we need to split the vector into two vectors first 
                // then borrow each node mutably
                //
                // it will be the last element of the first slice 
                // and first element of the second slice

                let (vec_slice_one, vec_slice_two) = 
                fluid_vec_in_loop.split_at_mut(index+1);

                let fluid_node_idx: &mut HeatTransferEntity = 
                vec_slice_one.last_mut().unwrap();

                let fluid_node_idx_plus_one: &mut HeatTransferEntity = 
                vec_slice_two.first_mut().unwrap();

                // after doing all the acrobatics to borrow two vectors, 
                // then we get the densities

                let fluid_node_idx_temperature: ThermodynamicTemperature = 
                HeatTransferEntity::temperature(
                    fluid_node_idx
                ).unwrap();

                let fluid_node_idx_plus_one_temperature: 
                ThermodynamicTemperature = 
                HeatTransferEntity::temperature(
                    fluid_node_idx_plus_one
                ).unwrap();


                let fluid_node_idx_density = density(
                    therminol,
                    fluid_node_idx_temperature,
                    atmospheric_pressure
                ).unwrap();

                let fluid_node_idx_plus_one_density = density(
                    therminol,
                    fluid_node_idx_plus_one_temperature,
                    atmospheric_pressure
                ).unwrap();

                // construct the advection interaction

                let mid_heater_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: fluid_node_idx_density,
                    fluid_density_heat_transfer_entity_2: fluid_node_idx_plus_one_density,
                };


                let mid_heater_advection_interaction = HeatTransferInteractionType::
                    Advection(mid_heater_advection_dataset);
                
                // link the nodes with advection

                link_heat_transfer_entity(fluid_node_idx, 
                    fluid_node_idx_plus_one, 
                    mid_heater_advection_interaction).unwrap();

            }

            // fifth, link fluid boundary nodes with the boundary 
            // conditions

            {
                // inlet link
                let therminol_inlet_temperature = 
                ThermodynamicTemperature::new::<degree_celsius>(80.0);

                let therminol_inlet_density = density(
                    therminol,
                    therminol_inlet_temperature,
                    atmospheric_pressure
                ).unwrap();

                // i need to borrow the first and last indexed vector 
                // mutably, so I split again

                let (vec_slice_one, vec_slice_two) = 
                fluid_vec_in_loop.split_at_mut(1);

                let mut inlet_node: &mut HeatTransferEntity = 
                vec_slice_one.first_mut().unwrap();

                let mut outlet_node: &mut HeatTransferEntity = 
                vec_slice_two.last_mut().unwrap();

                // we now need inlet and outlet node densities

                let inlet_node_density_vec = 
                HeatTransferEntity::density_vector( 
                    inlet_node.deref_mut()).unwrap();

                let inlet_node_density: MassDensity = 
                inlet_node_density_vec[0];

                let outlet_node_density_vec = 
                HeatTransferEntity::density_vector( 
                    outlet_node.deref_mut()).unwrap();

                let outlet_node_density: MassDensity = 
                outlet_node_density_vec[0];

                // construct the interaction objects to say we 
                // have an advection going on between the nodes and the 
                // BCs
                // then we can make the interactions work

                let inlet_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: therminol_inlet_density,
                    fluid_density_heat_transfer_entity_2: inlet_node_density,
                };

                let outlet_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: outlet_node_density,
                    // cv2 doesn't really matter here,
                    fluid_density_heat_transfer_entity_2: outlet_node_density,
                };


                let inlet_interaction = HeatTransferInteractionType::
                    Advection(inlet_advection_dataset);
                let outlet_interaction = HeatTransferInteractionType::
                    Advection(outlet_advection_dataset);


                // link the inlet and outlet with their respective BCs 
                //
                // (inlet bc) --- (inlet node) --- ... --- (outlet node) --- (outlet bc)

                link_heat_transfer_entity(&mut inlet_const_temp_in_loop, 
                    &mut inlet_node, 
                    inlet_interaction).unwrap();

                link_heat_transfer_entity(&mut outlet_node, 
                    &mut outlet_zero_heat_flux_in_loop, 
                    outlet_interaction).unwrap();

            }

            let node_connection_end_ms = 
            node_connection_start.elapsed().unwrap().as_millis();


            // I also want to see what the automatic timestepping 
            // is 


            // todo:
            {

                let mut temp_profile_data_vec: Vec<String> = vec![];
                // code block for recording temperature profiles
                // across node surfaces
                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                // next simulation time string 
                let elapsed_calc_time_seconds_string = 
                calculation_time_elapsed.elapsed().unwrap().as_secs().to_string();

                temp_profile_data_vec.push(current_time_string);
                temp_profile_data_vec.push(elapsed_calc_time_seconds_string);

                for outer_shell_ptr in steel_shell_outer_node_vec_in_loop.iter_mut(){

                    // get temperature first 
                    
                    let node_temp: ThermodynamicTemperature = 
                    HeatTransferEntity::temperature( 
                        outer_shell_ptr).unwrap();

                    // get it in degc 

                    let node_temp_deg_c: f64 = 
                    node_temp.get::<degree_celsius>();

                    // convert to string and push to vector 

                    let node_temp_c_string: String = 
                    node_temp_deg_c.to_string();

                    temp_profile_data_vec.push(node_temp_c_string);

                }

                // now write 
                temp_profile_wtr.write_record(&temp_profile_data_vec).unwrap();
                

            }
            
            let data_recording_time_start = SystemTime::now();


            
            let mut auto_calculated_timestep: Time = Time::new::<second>(100.0);
            // todo:
            {
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

                let st_11_node: &mut HeatTransferEntity = 
                &mut steel_shell_outer_node_vec_in_loop[st_11_index_number];

                // now i also want the therminol outlet temperature 

                let thermoinol_outlet_node: &mut HeatTransferEntity = 
                fluid_vec_in_loop.last_mut().unwrap();

                let therminol_outlet_temp_string = 
                HeatTransferEntity::temperature(
                    thermoinol_outlet_node).unwrap()
                    .get::<degree_celsius>().to_string();

                let shell_celsius_string = 
                HeatTransferEntity::temperature(
                    st_11_node).unwrap()
                    .get::<degree_celsius>().to_string();

                // drop the mutable references manually

                drop(thermoinol_outlet_node);
                drop(st_11_node);

                // then I want to get the max timestep
                //
                // 
                // get max timestep with 5 C max temp change 
                // this will usually ensure that max timestep change is 
                // dependent on Co and Fo, 
                //
                // I'll need to loop over all values to get the correct 
                // timestep


                // now loop over all nodes, get its temperature 

                for therminol_node in fluid_vec_in_loop.iter_mut() {


                    let local_timestep: Time = HeatTransferEntity::
                        get_max_timestep(therminol_node,
                            TemperatureInterval::new::<
                                uom::si::temperature_interval::kelvin>(5.0))
                        .unwrap();

                    if local_timestep < auto_calculated_timestep {
                        auto_calculated_timestep = local_timestep
                    }

                }

                for steel_inner_node in steel_shell_inner_node_vec_in_loop.iter_mut() {


                    let local_timestep: Time = HeatTransferEntity::
                        get_max_timestep(steel_inner_node,
                            TemperatureInterval::new::<
                                uom::si::temperature_interval::kelvin>(5.0))
                        .unwrap();

                    if local_timestep < auto_calculated_timestep {
                        auto_calculated_timestep = local_timestep
                    }

                }


                for steel_outer_node in steel_shell_outer_node_vec_in_loop.iter_mut() {


                    let local_timestep: Time = HeatTransferEntity::
                        get_max_timestep(steel_outer_node,
                            TemperatureInterval::new::<
                                uom::si::temperature_interval::kelvin>(5.0))
                        .unwrap();

                    if local_timestep < auto_calculated_timestep {
                        auto_calculated_timestep = local_timestep
                    }

                }

                
                let auto_calculated_timestep_string = 
                auto_calculated_timestep.get::<second>().to_string();


                // csv data writing
                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                let heater_power_kilowatt_string = 
                heater_power.get::<kilowatt>().to_string();

                wtr.write_record(&[current_time_string,
                    heater_power_kilowatt_string,
                    therminol_outlet_temp_string,
                    shell_celsius_string,
                    auto_calculated_timestep_string])
                    .unwrap();


            }
            let data_recording_time_end_ms = 
            data_recording_time_start.elapsed().unwrap().as_millis();

            let timestep_advance_start = 
            SystemTime::now();

            // advancing timestep, 
            // this part is extremely sluggish. However, it is paralellisable 
            // need to use Rayon here, otherwise, it will take forever
            //
            // or else some other parallelism here
            {
                // code block for advancing timestep over all control 
                // volumes, that is inner shell, outer shell and 
                // fluid volumes

                for therminol_node in fluid_vec_in_loop.iter_mut() {
                    HeatTransferEntity::advance_timestep(
                        therminol_node,
                        auto_calculated_timestep).unwrap();

                }

                for steel_inner_node in steel_shell_inner_node_vec_in_loop.iter_mut() {

                    HeatTransferEntity::advance_timestep(
                        steel_inner_node,
                        auto_calculated_timestep).unwrap();

                }


                for steel_outer_node in steel_shell_outer_node_vec_in_loop.iter_mut() {

                    HeatTransferEntity::advance_timestep(
                        steel_outer_node,
                        auto_calculated_timestep).unwrap();

                }
            }


            let timestep_advance_end_ms = 
            timestep_advance_start.elapsed().unwrap().as_millis();

            // write timestep diagnostics

            let total_time_ms = 
            arc_mutex_lock_elapsed_ms 
            + node_connection_end_ms 
            + data_recording_time_end_ms 
            + timestep_advance_end_ms;

            let mutex_lock_frac: f64 = arc_mutex_lock_elapsed_ms as f64 / 
            total_time_ms as f64;

            let node_connection_frac: f64 = node_connection_end_ms as f64 / 
            total_time_ms as f64;

            let data_record_frac: f64 = data_recording_time_end_ms as f64 / 
            total_time_ms as f64; 

            let timestep_advance_frac: f64 = timestep_advance_end_ms as f64 / 
            total_time_ms as f64;

            time_wtr.write_record(&[total_time_ms.to_string(),
                mutex_lock_frac.to_string(),
                node_connection_frac.to_string(),
                data_record_frac.to_string(),
                timestep_advance_frac.to_string()])
                .unwrap();
            

            // add the timestep
            current_time_simulation_time += auto_calculated_timestep;
        }
        // with csvs being written,
        // use cargo watch -x test --ignore '*.csv'
        wtr.flush().unwrap();
        time_wtr.flush().unwrap();
        temp_profile_wtr.flush().unwrap();
    };

    let calculation_thread = thread::spawn(calculation_loop);
    
    calculation_thread.join().unwrap();

    // done!
    return ();

}


/// This tests how reducing to one radial steel node reduces required 
/// time step
#[test]
pub fn ciet_heater_v_2_0_test_steady_state_functional_test_v_1_2_one_steel_node(){


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

    // construct the objects,
    // I'm going to use a function 

    fn construct_heated_section_fluid_nodes(therminol: Material,
        cross_sectional_area: Area,
        heated_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let mut fluid_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = heated_length/number_of_nodes as f64;

        for _index in 0..number_of_nodes {
            let therminol_node: HeatTransferEntity = 
            SingleCVNode::new_odd_shaped_pipe(
                node_length,
                cross_sectional_area,
                therminol,
                initial_temperature,
                pressure,
            ).unwrap();

            fluid_node_vec.push(therminol_node);
        }

        return fluid_node_vec;
    }

    let fluid_node_vec: Vec<HeatTransferEntity> 
    = construct_heated_section_fluid_nodes(
        therminol,
        flow_area,
        heated_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);
    
    // then construct two layers of steel shells
    
    fn construct_steel_shell_nodes(steel: Material,
        id: Length, 
        od: Length,
        heated_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let id: InnerDiameterThermalConduction = id.into();
        let od: OuterDiameterThermalConduction = od.into();

        let mut steel_shell_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = heated_length/number_of_nodes as f64;

        for _index in 0..number_of_nodes {

            let steel_shell_node = SingleCVNode::new_cylindrical_shell(
                node_length,
                id, od,
                steel, 
                initial_temperature,
                pressure,
            ).unwrap();

            steel_shell_node_vec.push(steel_shell_node);
        }

        return steel_shell_node_vec;
    }

    // inner layer of steel shell I will just assume 0.0392 m is the 
    // midway point
    // the inner node should be thicker anyway 


    let steel_shell_inner_node_vec: Vec<HeatTransferEntity> = 
    construct_steel_shell_nodes(
        steel,
        id, od, heated_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);


    // now, let me make mutex locks and Arc pointers
    // for the fluid nodes and steel nodes

    let fluid_node_vec_ptr = Arc::new(Mutex::new(
        fluid_node_vec
    ));

    let steel_shell_inner_node_vec_ptr = Arc::new(Mutex::new(
        steel_shell_inner_node_vec
    ));

    

    // need two boundary conditions 

    let inlet_const_temp = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ThermodynamicTemperature::new::<degree_celsius>(79.12)
        ));

    let outlet_zero_heat_flux = HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatAddition(Power::new::<watt>(0.0))
    );

    let ambient_temperature_bc = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ambient_air_temp
        ));

    let inlet_const_temp_ptr = Arc::new(Mutex::new(
        inlet_const_temp
    ));

    let outlet_zero_heat_flux_ptr = Arc::new(Mutex::new(
        outlet_zero_heat_flux
    ));

    let ambient_air_temp_bc_ptr = Arc::new(Mutex::new(
        ambient_temperature_bc
    ));

    // the two types of HeatTransferInteractionType are 
    // advection and convection resistance
    //
    // 2007 square_centimeter
    // and 607 watt_per_square_meter_kelvin



    // timestep settings


    let max_time: Time = Time::new::<second>(200.0);
    let max_time_ptr = Arc::new(max_time);

    let calculation_time_elapsed = SystemTime::now();

    // this is the calculation loop
    let calculation_loop = move || {


        // csv writer, for post processing 

        // now for heater power, can swap between 

        let heater_power = heater_steady_state_power;


        let mut wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_one_steel_node.csv")
            .unwrap();

        wtr.write_record(&["time_seconds",
            "heater_power_kilowatts",
            "therminol_temperature_celsius",
            "shell_temperature_celsius",
            "auto_timestep_calculated_seconds",])
            .unwrap();

        let mut time_wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_time_one_steel_node.csv")
            .unwrap();

        time_wtr.write_record(&["loop_calculation_time_ms",
            "mutex_lock_frac",
            "node_connection_frac",
            "data_record_frac",
            "timestep_advance_frac",])
            .unwrap();

        // now i want a writer for temperature profile over n nodes 
        // it needs a simulation time, computation time elapsed, 
        // and a temperature for the outer surface temperature node 
        // for all nodes 

        let mut temp_profile_wtr = Writer::from_path("ciet_heater_v_2_0_steady_state_temp_profile_one_steel_node.csv")
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


        let mut current_time_simulation_time = Time::new::<second>(0.0);

        let max_time_ptr_in_loop = max_time_ptr;
        // we are sampling at about 10 Hz
        // so the nyquist frequency is about 5 Hz 
        // this is because the highest frequency is about 3.66 Hz
        
        while current_time_simulation_time <= *max_time_ptr_in_loop {

            // timer for timekeeping purposes 

            let arc_mutex_lock_start = SystemTime::now();

            // calculation steps

            let mut inlet_const_temp_in_loop = 
            inlet_const_temp_ptr.lock().unwrap();
            let mut outlet_zero_heat_flux_in_loop = 
            outlet_zero_heat_flux_ptr.lock().unwrap();
            
            let mut ambient_air_temp_bc_in_loop = 
            ambient_air_temp_bc_ptr.lock().unwrap();

            let mut fluid_vec_in_loop = 
            fluid_node_vec_ptr.lock().unwrap();

            let mut steel_shell_inner_node_vec_in_loop = 
            steel_shell_inner_node_vec_ptr.lock().unwrap();


            let arc_mutex_lock_end = SystemTime::now();

            let arc_mutex_lock_elapsed_ms = 
            arc_mutex_lock_start.elapsed().unwrap().as_millis()
            - arc_mutex_lock_end.elapsed().unwrap().as_millis();

            // we need to conenct a few things 
            //
            // to simplify, we also ignore axial conduction 
            // in all materials
            // 
            // in the radial direction:
            //
            // (ambient temp)
            // |
            // | --  convection/conduction
            // |
            // (steel shell nodes) 
            // | 
            // | -- conduction/convection
            // | 
            // (fluid nodes) 
            //
            // 
            // in the axial direction: 
            //
            // (inlet) --- (fluid nodes) --- (outlet)

            // calculate convection interactions
            // first, we need to connect the 
            //

            let node_connection_start = SystemTime::now();


            // radial heat transfer interactions
            // this one is fluid to inner steel cylinder
            for (index, fluid_node_heat_trf_entity_ptr) in 
                fluid_vec_in_loop.iter_mut().enumerate(){

                    // conduction material, properties

                    // radial thickness needs to be half of the 
                    // shell thickness, because it links to the shell 
                    // center
                    let radial_thickness: Length = 
                    (od - id) *0.5;

                    let radial_thickness: RadialCylindricalThicknessThermalConduction
                    = radial_thickness.into();

                    let steel_inner_cylindrical_node_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature(
                        &mut steel_shell_inner_node_vec_in_loop[index]).unwrap();

                    let pressure = atmospheric_pressure;

                    // now need to get heat transfer coeff

                    let therminol_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature(
                        fluid_node_heat_trf_entity_ptr).unwrap();

                    let heat_trf_coeff: HeatTransfer = 
                    heat_transfer_coefficient_ciet_v_2_0(
                        therminol_mass_flowrate,
                        therminol_temp,
                        pressure,
                    );

                    let inner_diameter: InnerDiameterThermalConduction = 
                    id.clone().into();

                    let node_length: Length = 
                    heated_length/(number_of_nodes as f64);

                    let node_length: CylinderLengthThermalConduction = 
                    node_length.into();

                    // construct the interaction 

                    let interaction: HeatTransferInteractionType = 
                    HeatTransferInteractionType::CylindricalConductionConvectionLiquidInside
                        ((steel,radial_thickness,
                            steel_inner_cylindrical_node_temp,
                            pressure),
                            (heat_trf_coeff,
                                inner_diameter,
                                node_length));

                    // link the entities,
                    // this is the fluid to the inner shell
                    link_heat_transfer_entity(fluid_node_heat_trf_entity_ptr, 
                        &mut steel_shell_inner_node_vec_in_loop[index], 
                        interaction).unwrap();

                }

            // third, link outer shell to ambient temperature

            for (_index,outer_shell_ptr) in 
                steel_shell_inner_node_vec_in_loop.iter_mut().enumerate(){

                    // conduction material, properties

                    // radial thickness needs to be half of the 
                    // shell thickness, because it links to the shell 
                    // center
                    //
                    // (outer shell) ----- (outer shell surface) ---- (fluid)
                    // 
                    //          R_{shell} (half length)    R_{conv}
                    //
                    //
                    //  for convection, h = 20 W/(m^2 K)
                    //
                    let radial_thickness: Length = 
                    (od - id) *0.5;

                    let radial_thickness: RadialCylindricalThicknessThermalConduction
                    = radial_thickness.into();

                    let steel_inner_cylindrical_node_temp: ThermodynamicTemperature 
                    = HeatTransferEntity::temperature( 
                        outer_shell_ptr).unwrap();

                    let pressure = atmospheric_pressure;

                    // now need to get heat transfer coeff
                    //
                    // 20 W/(m^2 K)

                    let heat_trf_coeff: HeatTransfer = 
                    HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);


                    let outer_diameter: OuterDiameterThermalConduction = 
                    od.clone().into();

                    let node_length: Length = 
                    heated_length/(number_of_nodes as f64);

                    let node_length: CylinderLengthThermalConduction = 
                    node_length.into();

                    // construct the interaction 

                    let interaction: HeatTransferInteractionType = 
                    HeatTransferInteractionType::CylindricalConductionConvectionLiquidOutside
                        ((steel,radial_thickness,
                            steel_inner_cylindrical_node_temp,
                            pressure),
                            (heat_trf_coeff,
                                outer_diameter,
                                node_length));

                    // link the entities,
                    // this is the fluid to the inner shell
                    link_heat_transfer_entity(&mut ambient_air_temp_bc_in_loop, 
                        outer_shell_ptr, 
                        interaction).unwrap();

                    


                }
            
            // fourth, link adjacent fluid nodes axially with advection 

            for index in 0..fluid_vec_in_loop.len()-1 {

                // (heater node i) ----- (heater node 2) ----- .... 
                //
                //
                //
                // to borrow two elements mutably is tricky, 
                // so we need to split the vector into two vectors first 
                // then borrow each node mutably
                //
                // it will be the last element of the first slice 
                // and first element of the second slice

                let (vec_slice_one, vec_slice_two) = 
                fluid_vec_in_loop.split_at_mut(index+1);

                let fluid_node_idx: &mut HeatTransferEntity = 
                vec_slice_one.last_mut().unwrap();

                let fluid_node_idx_plus_one: &mut HeatTransferEntity = 
                vec_slice_two.first_mut().unwrap();

                // after doing all the acrobatics to borrow two vectors, 
                // then we get the densities

                let fluid_node_idx_temperature: ThermodynamicTemperature = 
                HeatTransferEntity::temperature(
                    fluid_node_idx
                ).unwrap();

                let fluid_node_idx_plus_one_temperature: 
                ThermodynamicTemperature = 
                HeatTransferEntity::temperature(
                    fluid_node_idx_plus_one
                ).unwrap();


                let fluid_node_idx_density = density(
                    therminol,
                    fluid_node_idx_temperature,
                    atmospheric_pressure
                ).unwrap();

                let fluid_node_idx_plus_one_density = density(
                    therminol,
                    fluid_node_idx_plus_one_temperature,
                    atmospheric_pressure
                ).unwrap();

                // construct the advection interaction

                let mid_heater_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: fluid_node_idx_density,
                    fluid_density_heat_transfer_entity_2: fluid_node_idx_plus_one_density,
                };


                let mid_heater_advection_interaction = HeatTransferInteractionType::
                    Advection(mid_heater_advection_dataset);
                
                // link the nodes with advection

                link_heat_transfer_entity(fluid_node_idx, 
                    fluid_node_idx_plus_one, 
                    mid_heater_advection_interaction).unwrap();

            }

            // fifth, link fluid boundary nodes with the boundary 
            // conditions

            {
                // inlet link
                let therminol_inlet_temperature = 
                ThermodynamicTemperature::new::<degree_celsius>(80.0);

                let therminol_inlet_density = density(
                    therminol,
                    therminol_inlet_temperature,
                    atmospheric_pressure
                ).unwrap();

                // i need to borrow the first and last indexed vector 
                // mutably, so I split again

                let (vec_slice_one, vec_slice_two) = 
                fluid_vec_in_loop.split_at_mut(1);

                let mut inlet_node: &mut HeatTransferEntity = 
                vec_slice_one.first_mut().unwrap();

                let mut outlet_node: &mut HeatTransferEntity = 
                vec_slice_two.last_mut().unwrap();

                // we now need inlet and outlet node densities

                let inlet_node_density_vec = 
                HeatTransferEntity::density_vector( 
                    inlet_node.deref_mut()).unwrap();

                let inlet_node_density: MassDensity = 
                inlet_node_density_vec[0];

                let outlet_node_density_vec = 
                HeatTransferEntity::density_vector( 
                    outlet_node.deref_mut()).unwrap();

                let outlet_node_density: MassDensity = 
                outlet_node_density_vec[0];

                // construct the interaction objects to say we 
                // have an advection going on between the nodes and the 
                // BCs
                // then we can make the interactions work

                let inlet_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: therminol_inlet_density,
                    fluid_density_heat_transfer_entity_2: inlet_node_density,
                };

                let outlet_advection_dataset = DataAdvection {
                    mass_flowrate: therminol_mass_flowrate,
                    fluid_density_heat_transfer_entity_1: outlet_node_density,
                    // cv2 doesn't really matter here,
                    fluid_density_heat_transfer_entity_2: outlet_node_density,
                };


                let inlet_interaction = HeatTransferInteractionType::
                    Advection(inlet_advection_dataset);
                let outlet_interaction = HeatTransferInteractionType::
                    Advection(outlet_advection_dataset);


                // link the inlet and outlet with their respective BCs 
                //
                // (inlet bc) --- (inlet node) --- ... --- (outlet node) --- (outlet bc)

                link_heat_transfer_entity(&mut inlet_const_temp_in_loop, 
                    &mut inlet_node, 
                    inlet_interaction).unwrap();

                link_heat_transfer_entity(&mut outlet_node, 
                    &mut outlet_zero_heat_flux_in_loop, 
                    outlet_interaction).unwrap();

            }

            let node_connection_end_ms = 
            node_connection_start.elapsed().unwrap().as_millis();


            // I also want to see what the automatic timestepping 
            // is 


            // todo:
            {

                let mut temp_profile_data_vec: Vec<String> = vec![];
                // code block for recording temperature profiles
                // across node surfaces
                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                // next simulation time string 
                let elapsed_calc_time_seconds_string = 
                calculation_time_elapsed.elapsed().unwrap().as_secs().to_string();

                temp_profile_data_vec.push(current_time_string);
                temp_profile_data_vec.push(elapsed_calc_time_seconds_string);

                for inner_shell_ptr in steel_shell_inner_node_vec_in_loop.iter_mut(){

                    // get temperature first 
                    
                    let node_temp: ThermodynamicTemperature = 
                    HeatTransferEntity::temperature( 
                        inner_shell_ptr).unwrap();

                    // get it in degc 

                    let node_temp_deg_c: f64 = 
                    node_temp.get::<degree_celsius>();

                    // convert to string and push to vector 

                    let node_temp_c_string: String = 
                    node_temp_deg_c.to_string();

                    temp_profile_data_vec.push(node_temp_c_string);

                }

                // now write 
                temp_profile_wtr.write_record(&temp_profile_data_vec).unwrap();
                

            }
            
            let data_recording_time_start = SystemTime::now();


            
            let mut auto_calculated_timestep: Time = Time::new::<second>(100.0);
            // todo:
            {
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

                let st_11_node: &mut HeatTransferEntity = 
                &mut steel_shell_inner_node_vec_in_loop[st_11_index_number];

                // now i also want the therminol outlet temperature 

                let thermoinol_outlet_node: &mut HeatTransferEntity = 
                fluid_vec_in_loop.last_mut().unwrap();

                let therminol_outlet_temp_string = 
                HeatTransferEntity::temperature(
                    thermoinol_outlet_node).unwrap()
                    .get::<degree_celsius>().to_string();

                let shell_celsius_string = 
                HeatTransferEntity::temperature(
                    st_11_node).unwrap()
                    .get::<degree_celsius>().to_string();

                // drop the mutable references manually

                drop(thermoinol_outlet_node);
                drop(st_11_node);

                // then I want to get the max timestep
                //
                // 
                // get max timestep with 5 C max temp change 
                // this will usually ensure that max timestep change is 
                // dependent on Co and Fo, 
                //
                // I'll need to loop over all values to get the correct 
                // timestep


                // now loop over all nodes, get its temperature 

                for therminol_node in fluid_vec_in_loop.iter_mut() {


                    let local_timestep: Time = HeatTransferEntity::
                        get_max_timestep(therminol_node,
                            TemperatureInterval::new::<
                                uom::si::temperature_interval::kelvin>(5.0))
                        .unwrap();

                    if local_timestep < auto_calculated_timestep {
                        auto_calculated_timestep = local_timestep
                    }

                }

                for steel_inner_node in steel_shell_inner_node_vec_in_loop.iter_mut() {


                    let local_timestep: Time = HeatTransferEntity::
                        get_max_timestep(steel_inner_node,
                            TemperatureInterval::new::<
                                uom::si::temperature_interval::kelvin>(5.0))
                        .unwrap();

                    if local_timestep < auto_calculated_timestep {
                        auto_calculated_timestep = local_timestep
                    }

                }



                
                let auto_calculated_timestep_string = 
                auto_calculated_timestep.get::<second>().to_string();


                // csv data writing
                let current_time_string = 
                current_time_simulation_time.get::<second>().to_string();

                let heater_power_kilowatt_string = 
                heater_power.get::<kilowatt>().to_string();

                wtr.write_record(&[current_time_string,
                    heater_power_kilowatt_string,
                    therminol_outlet_temp_string,
                    shell_celsius_string,
                    auto_calculated_timestep_string])
                    .unwrap();


            }
            let data_recording_time_end_ms = 
            data_recording_time_start.elapsed().unwrap().as_millis();

            let timestep_advance_start = 
            SystemTime::now();

            // advancing timestep, 
            // this part is extremely sluggish. However, it is paralellisable 
            // need to use Rayon here, otherwise, it will take forever
            //
            // or else some other parallelism here
            {
                // code block for advancing timestep over all control 
                // volumes, that is inner shell, outer shell and 
                // fluid volumes

                for therminol_node in fluid_vec_in_loop.iter_mut() {
                    HeatTransferEntity::advance_timestep(
                        therminol_node,
                        auto_calculated_timestep).unwrap();

                }

                for steel_inner_node in steel_shell_inner_node_vec_in_loop.iter_mut() {

                    HeatTransferEntity::advance_timestep(
                        steel_inner_node,
                        auto_calculated_timestep).unwrap();

                }


            }


            let timestep_advance_end_ms = 
            timestep_advance_start.elapsed().unwrap().as_millis();

            // write timestep diagnostics

            let total_time_ms = 
            arc_mutex_lock_elapsed_ms 
            + node_connection_end_ms 
            + data_recording_time_end_ms 
            + timestep_advance_end_ms;

            let mutex_lock_frac: f64 = arc_mutex_lock_elapsed_ms as f64 / 
            total_time_ms as f64;

            let node_connection_frac: f64 = node_connection_end_ms as f64 / 
            total_time_ms as f64;

            let data_record_frac: f64 = data_recording_time_end_ms as f64 / 
            total_time_ms as f64; 

            let timestep_advance_frac: f64 = timestep_advance_end_ms as f64 / 
            total_time_ms as f64;

            time_wtr.write_record(&[total_time_ms.to_string(),
                mutex_lock_frac.to_string(),
                node_connection_frac.to_string(),
                data_record_frac.to_string(),
                timestep_advance_frac.to_string()])
                .unwrap();
            

            // add the timestep
            current_time_simulation_time += auto_calculated_timestep;
        }
        // with csvs being written,
        // use cargo watch -x test --ignore '*.csv'
        wtr.flush().unwrap();
        time_wtr.flush().unwrap();
        temp_profile_wtr.flush().unwrap();
    };

    let calculation_thread = thread::spawn(calculation_loop);
    
    calculation_thread.join().unwrap();

    // done!
    return ();


}

fn heat_transfer_coefficient_ciet_v_2_0(mass_flowrate: MassRate,
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


    let reynolds: Ratio = ciet_heater_v_2_0_reynolds_nunber(
        mass_flowrate, mu);

    let prandtl: Ratio = liquid_prandtl(
        therminol,
        therminol_temperature,
        pressure).unwrap();

    let nusselt: Ratio = ciet_heater_v_2_0_nusselt_number(
        reynolds,
        prandtl);

    let hydraulic_diameter = Length::new::<meter>(0.01467);

    let heat_transfer_coeff: HeatTransfer = 
    nusselt * k / hydraulic_diameter;

    heat_transfer_coeff

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

#[inline]
fn ciet_heater_v_2_0_reynolds_nunber(mass_flowrate: MassRate,
    mu: DynamicViscosity) -> Ratio {

    // Re = m* D_H/ A_{XS}/mu
    let hydraulic_diameter = Length::new::<meter>(0.01467);
    let flow_area = Area::new::<square_meter>(0.00105);

    mass_flowrate*hydraulic_diameter/mu/flow_area
}
