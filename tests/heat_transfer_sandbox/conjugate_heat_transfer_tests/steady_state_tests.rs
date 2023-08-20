use std::ops::DerefMut;
use std::sync::{Arc, Mutex};
use std::thread;

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
use uom::si::length::{centimeter, meter};
use uom::si::mass_rate::kilogram_per_second;
use uom::si::power::{watt, kilowatt};
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::time::second;

use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations
::heat_transfer_entities::BCType::*;

/// In this test, we have a nodalised representation of the 
/// this is 8 nodes in the axial direction and 2 nodes for metal 
/// in the radial direction
///
/// This is for heater v2.0
///
///
#[test]
//#[ignore = "takes about 20min, only use for data collection"]
pub fn ciet_heater_v_1_0_test_steady_state(){


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
    let heated_length = Length::new::<meter>(1.676);
    let total_length = Length::new::<meter>(1.983333);
    let initial_temperature = ThermodynamicTemperature::new::
        <degree_celsius>(80.0);
    let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

    let flow_area = Area::new::<square_meter>(0.00105);
    let number_of_nodes: usize = 8;
    let ambient_air_temp = ThermodynamicTemperature::new::<
        degree_celsius>(21.67);

    
    

    let therminol_mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);

    // construct the objects,
    // I'm going to use a function 

    fn construct_heated_section_fluid_nodes(therminol: Material,
        cross_sectional_area: Area,
        total_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let mut fluid_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = total_length/number_of_nodes as f64;

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

    let mut fluid_node_vec: Vec<HeatTransferEntity> 
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
        total_length: Length,
        initial_temperature: ThermodynamicTemperature,
        pressure: Pressure,
        number_of_nodes: usize,) -> Vec<HeatTransferEntity>{

        // I'm going to make a vector of mutable heat transfer 
        // entities

        let id: InnerDiameterThermalConduction = id.into();
        let od: OuterDiameterThermalConduction = od.into();

        let mut steel_shell_node_vec: Vec<HeatTransferEntity> = vec![];

        // now let's get individual length of each node 

        let node_length: Length = total_length/number_of_nodes as f64;

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

    let mut steel_shell_inner_node_vec: Vec<HeatTransferEntity> = 
    construct_steel_shell_nodes(
        steel,
        id, midway_point_steel_shell, total_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);

    let mut steel_shell_outer_node_vec: Vec<HeatTransferEntity> = 
    construct_steel_shell_nodes(
        steel,
        midway_point_steel_shell, od, total_length,
        initial_temperature,
        atmospheric_pressure,
        number_of_nodes);

    // create mutex locks for them 

    let fluid_node_vec_ptr = Arc::new(Mutex::new(
        fluid_node_vec
    ));

    let steel_shell_inner_node_vec_ptr = Arc::new(Mutex::new(
        steel_shell_inner_node_vec
    ));

    let steel_shell_outer_node_vec_ptr = Arc::new(Mutex::new(
        steel_shell_outer_node_vec
    ));
    

    let therminol_cylinder: HeatTransferEntity = 
    SingleCVNode::new_cylindrical_shell(
        total_length,
        inner_tube_od.into(),
        id.into(),
        therminol,
        initial_temperature,
        atmospheric_pressure
    ).unwrap();

    let steel_shell: HeatTransferEntity = 
    SingleCVNode::new_cylindrical_shell(
        heated_length,
        id.into(),
        od.into(),
        steel,
        initial_temperature,
        atmospheric_pressure
    ).unwrap();

    // now, let me make mutex locks and Arc pointers

    let therminol_cylinder_ptr = Arc::new(
        Mutex::new(
            therminol_cylinder
        ));

    let steel_shell_ptr = Arc::new(
        Mutex::new(
            steel_shell
        ));



    // need two boundary conditions 

    let inlet_const_temp = HeatTransferEntity::BoundaryConditions(
        UserSpecifiedTemperature(
            ThermodynamicTemperature::new::<degree_celsius>(80.0)
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


    let max_time: Time = Time::new::<second>(2.0);
    let max_time_ptr = Arc::new(max_time);

    // this is the calculation loop
    let calculation_loop = move || {
        // csv writer, for post processing 


        let mut wtr = Writer::from_path("one_dimension_ciet_cht.csv")
            .unwrap();

        wtr.write_record(&["time_seconds",
            "heater_power_kilowatts",
            "therminol_temperature_celsius",
            "shell_temperature_celsius",
            "auto_timestep_calculated_seconds",])
            .unwrap();


        let mut current_time_simulation_time = Time::new::<second>(0.0);

        let max_time_ptr_in_loop = max_time_ptr;
        // we are sampling at about 10 Hz
        // so the nyquist frequency is about 5 Hz 
        // this is because the highest frequency is about 3.66 Hz
        let timestep_value = Time::new::<second>(0.1);
        
        while current_time_simulation_time <= *max_time_ptr_in_loop {

            // calculation steps

            let mut therminol_cylinder_in_loop = 
            therminol_cylinder_ptr.lock().unwrap();
            let mut steel_shell_in_loop = 
            steel_shell_ptr.lock().unwrap();
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
                    total_length/(number_of_nodes as f64);

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
                    total_length/(number_of_nodes as f64);

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

                    // now for heater power 


                    let heater_power = mfbs_power_signal_logspace_custom(
                        current_time_simulation_time);

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
                    total_length/(number_of_nodes as f64);

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



            // advection bc, so at boundary condition, therminol flows in at 
            // 0.18 kg/s at 80C


            let therminol_inlet_temperature = 
            ThermodynamicTemperature::new::<degree_celsius>(80.0);



            let therminol_inlet_density = density(
                therminol,
                therminol_inlet_temperature,
                atmospheric_pressure
            ).unwrap();

            // i can also calculate the densities of each cv

            let therminol_cv_density_vec = 
            HeatTransferEntity::density_vector( 
                therminol_cylinder_in_loop.deref_mut()).unwrap();

            let heater_fluid_cv_density: MassDensity = 
            therminol_cv_density_vec[0];

            // now crate the dataset 
            //
            // the diagram is like so 
            //
            // (inlet) -------------> cv --------------> (outlet)
            //
            //          inlet_interaction   outlet_interaction 
            //
            // the arguments are placed in order of left to right:
            //
            // inlet_advection_interaction(inlet,cv)
            // outlet_advection_interaction(cv,outlet)

            let inlet_advection_dataset = DataAdvection {
                mass_flowrate: therminol_mass_flowrate,
                fluid_density_heat_transfer_entity_1: therminol_inlet_density,
                fluid_density_heat_transfer_entity_2: heater_fluid_cv_density,
            };

            let outlet_advection_dataset = DataAdvection {
                mass_flowrate: therminol_mass_flowrate,
                fluid_density_heat_transfer_entity_1: heater_fluid_cv_density,
                // cv2 doesn't really matter here,
                fluid_density_heat_transfer_entity_2: therminol_inlet_density,
            };

            let inlet_interaction = HeatTransferInteractionType::
                Advection(inlet_advection_dataset);
            let outlet_interaction = HeatTransferInteractionType::
                Advection(outlet_advection_dataset);

            // now let's link the cv 

            // Electrical Heat
            // -------------> (solid shell) 
            //                     | 
            //                     | 
            //                     |  (thermal resistance)
            //                     | 
            //
            // --------> (well mixed fluid volume) ----------->
            //                 T_fluid
            // Fluid in                                Fluid out 
            // T_in                                    T_fluid

            // (1) inlet to fluid
            link_heat_transfer_entity(&mut inlet_const_temp_in_loop, 
                &mut therminol_cylinder_in_loop, 
                inlet_interaction).unwrap();
            // (2) fluid to outlet
            link_heat_transfer_entity(&mut therminol_cylinder_in_loop, 
                &mut outlet_zero_heat_flux_in_loop, 
                outlet_interaction).unwrap();

            // (3) therminol cylinder and steel shell cv
            link_heat_transfer_entity(&mut therminol_cylinder_in_loop, 
                &mut steel_shell_in_loop, 
                convection_resistance).unwrap();

            // (4) electrical heat to solid shell 
            // (todo)
            // will need to use the mfbs signal 

            let _heater_power = mfbs_poresky_2017_power_signal(
                current_time_simulation_time);

            let heater_power = mfbs_power_signal_logspace_custom(
                current_time_simulation_time);

            let mut electrical_heat_bc: HeatTransferEntity = 
            BCType::new_const_heat_addition(heater_power);

            let heat_addition_interaction = 
            HeatTransferInteractionType::UserSpecifiedHeatAddition;

            link_heat_transfer_entity(
                &mut steel_shell_in_loop,
                &mut electrical_heat_bc,
                heat_addition_interaction,
            ).unwrap();


            // I also want to see what the automatic timestepping 
            // is 

            let auto_calculated_timestep = HeatTransferEntity::
                get_max_timestep(&mut therminol_cylinder_in_loop,
                TemperatureInterval::new::<
                        uom::si::temperature_interval::kelvin>(20.0))
                .unwrap();

            // csv data writing
            let current_time_string = 
            current_time_simulation_time.get::<second>().to_string();

            let heater_power_kilowatt_string = 
            heater_power.get::<kilowatt>().to_string();

            let therminol_celsius_string = 
            HeatTransferEntity::temperature(
                &mut therminol_cylinder_in_loop).unwrap()
                .get::<degree_celsius>().to_string();

            let shell_celsius_string = 
            HeatTransferEntity::temperature(
                &mut steel_shell_in_loop).unwrap()
                .get::<degree_celsius>().to_string();


            let auto_calculated_timestep_string = 
            auto_calculated_timestep.get::<second>().to_string();



            wtr.write_record(&[current_time_string,
                heater_power_kilowatt_string,
                therminol_celsius_string,
                shell_celsius_string,
                auto_calculated_timestep_string])
                .unwrap();

            // advance timestep for steel shell and therminol cylinder
            HeatTransferEntity::advance_timestep(
                &mut therminol_cylinder_in_loop,
                timestep_value).unwrap();

            HeatTransferEntity::advance_timestep(
                &mut steel_shell_in_loop,
                timestep_value).unwrap();


            // add the timestep
            current_time_simulation_time += timestep_value;
        }
        // with csvs being written,
        // use cargo watch -x test --ignore '*.csv'
        wtr.flush().unwrap();
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
    // average it out (divide by 5)
    power_analog_signal *= 0.2;

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

fn analog_poresky_2017_power_signal(simulation_time: Time) -> Power {

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
    // average it out (divide by 5)
    0.2 * power_analog_signal

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

#[inline]
fn ciet_heater_v_2_0_reynolds_nunber(mass_flowrate: MassRate,
    mu: DynamicViscosity) -> Ratio {

    // Re = m* D_H/ A_{XS}/mu
    let hydraulic_diameter = Length::new::<meter>(0.01467);
    let flow_area = Area::new::<square_meter>(0.00105);

    mass_flowrate*hydraulic_diameter/mu/flow_area
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
