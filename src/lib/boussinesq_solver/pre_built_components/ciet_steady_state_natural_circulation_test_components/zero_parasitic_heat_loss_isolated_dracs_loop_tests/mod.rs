/// In the original SAM publication
///
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// I found it hard to distinguish what TCHX temperatures case A,
/// B and C were.
///
/// But there was another publication which shows which is test group 
/// corresponds to which temperature:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 3 also provides the data 
/// 
///
pub mod case_a;


/// In the original SAM publication
///
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// I found it hard to distinguish what TCHX temperatures case A,
/// B and C were.
///
/// But there was another publication which shows which is test group 
/// corresponds to which temperature:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 3 also provides the data 
/// 
pub mod case_b;

/// In the original SAM publication
///
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// I found it hard to distinguish what TCHX temperatures case A,
/// B and C were.
///
/// But there was another publication which shows which is test group 
/// corresponds to which temperature:
///
/// Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
/// SAM using natural-circulation experimental data from the compact 
/// integral effects test (CIET) facility. 
/// Nuclear Engineering and Design, 377, 111144.
///
/// According to table 2,
///
/// Case A has 7 tests and TCHX out temperature of 46 C
/// Case B has 9 tests and TCHX out temperature of 35 C
/// Case C has 9 tests and TCHX out temperature of 40 C
///
/// Table 3 also provides the data 
/// 
///

pub mod case_c;


/// debugging tests for thermal hydraulics and fluid mechanics 
/// functions to make natural circulation 
/// testing easier 
pub mod debugging_thermal_hydraulics;


/// debugging tests for PID controller
/// functions to make natural circulation 
/// testing easier 
pub mod debugging_pid_controller;

/// miscellaneous debugging tests 
/// for other bugs I happened to find
pub mod misc_debugging;


/// This next set of tests shows explicitly what we need to do in 
/// the fluid component collection in order to get natural circulation
///
/// prototype test two, 
///
/// found that I can't use closures woops
/// I'll just use simple functions instead 
/// quite lehcheh (troublesome) to type, but ok lah
///
/// probably need to check thermal resistance at the CTAH or NDHX 
/// I think the nusselt number correlation is a bottleneck because 
/// it is meant for a pipe, rather than a NDHX
///
/// In De Wet's work, I believe there was no fluid thermal resistance 
/// assumed in the heat exchanger (high nusselt number in other words)
///
#[test]
pub fn dracs_natural_circ_thermal_hydraulics_pid_test_prototype_2(){

    use uom::si::{frequency::hertz, ratio::ratio, time::millisecond};
    use uom::si::thermodynamic_temperature::kelvin;
    use uom::si::f64::*;

    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::dracs_loop_components::*;
    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_collection::FluidComponentCollection;
    // let's construct the branches with test pressures and obtain 
    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_collection::FluidComponentCollectionMethods;
    use uom::ConstZero;

    use uom::si::thermodynamic_temperature::degree_celsius;
    use crate::boussinesq_solver::
        array_control_vol_and_fluid_component_collections::
        fluid_component_collection::
        fluid_component_super_collection::FluidComponentSuperCollection;

    use crate::boussinesq_solver::pre_built_components::
        insulated_pipes_and_fluid_components::InsulatedFluidComponent;
    use crate::boussinesq_solver::pre_built_components::
        non_insulated_fluid_components::NonInsulatedFluidComponent;

    use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
    use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::
        heat_transfer_interaction_enums::HeatTransferInteractionType;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use uom::si::power::watt;
    use uom::si::time::second;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;


    use chem_eng_real_time_process_control_simulator::alpha_nightly::transfer_fn_wrapper_and_enums::TransferFnTraits;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::ProportionalController;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::AnalogController;
    // setup 
    let initial_temperature = ThermodynamicTemperature::
        new::<degree_celsius>(40.0);
    
    let timestep = Time::new::<second>(0.5);
    let heat_rate_through_dhx = Power::new::<watt>(931.8);
    let mut tchx_heat_transfer_coeff: HeatTransfer;

    let reference_tchx_htc = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(25.0);
    let average_temperature_for_density_calcs = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);
    // let's calculate 2000 seconds of simulated time 
    // it takes about that long for the temperature to settle down

    let mut current_simulation_time = Time::ZERO;
    let max_simulation_time = Time::new::<second>(2400.0);

    // PID controller settings
    let controller_gain = Ratio::new::<ratio>(1.75);
    let integral_time: Time = controller_gain / Frequency::new::<hertz>(1.0);
    let derivative_time: Time = Time::new::<second>(1.0);
    // derivative time ratio
    let alpha: Ratio = Ratio::new::<ratio>(1.0);

    let mut pid_controller: AnalogController = 
    AnalogController::new_filtered_pid_controller(controller_gain,
        integral_time,
        derivative_time,
        alpha).unwrap();

    // we also have a measurement delay of 0.0001 s 
    // or 0.1 ms
    let measurement_delay = Time::new::<millisecond>(0.1);

    let mut measurement_delay_block: AnalogController = 
    ProportionalController::new(Ratio::new::<ratio>(1.0)).unwrap().into();

    measurement_delay_block.set_dead_time(measurement_delay);

    // set point is 319 kelvin
    let tchx_outlet_temperature_set_point = 
        ThermodynamicTemperature::new::<degree_celsius>(46.0);

    // hot branch or (mostly) hot leg
    let mut pipe_34 = new_pipe_34(initial_temperature);
    let mut pipe_33 = new_pipe_33(initial_temperature);
    let mut pipe_32 = new_pipe_32(initial_temperature);
    let mut pipe_31a = new_pipe_31a(initial_temperature);
    let mut static_mixer_61_label_31 = new_static_mixer_61_label_31(initial_temperature);
    let mut dhx_tube_side_30b = new_dhx_tube_side_30b(initial_temperature);
    let mut dhx_tube_side_heat_exchanger_30 = new_isolated_dhx_tube_side_30(initial_temperature);
    let mut dhx_tube_side_30a = new_dhx_tube_side_30a(initial_temperature);

    // cold branch or (mostly) cold leg
    let mut tchx_35a = new_ndhx_tchx_horizontal_35a(initial_temperature);
    let mut tchx_35b = new_ndhx_tchx_vertical_35b(initial_temperature);
    let mut static_mixer_60_label_36 = new_static_mixer_60_label_36(initial_temperature);
    let mut pipe_36a = new_pipe_36a(initial_temperature);
    let mut pipe_37 = new_pipe_37(initial_temperature);
    let mut flowmeter_60_37a = new_flowmeter_60_37a(initial_temperature);
    let mut pipe_38 = new_pipe_38(initial_temperature);
    let mut pipe_39 = new_pipe_39(initial_temperature);

    // fluid mechanics bit 
    fn get_dracs_flowrate(dracs_branches: &FluidComponentSuperCollection) -> 
        MassRate {
            let pressure_change_across_each_branch = 
                dracs_branches.get_pressure_change(MassRate::ZERO);

            let mass_flowrate_across_each_branch: Vec<MassRate> = 
                dracs_branches.
                get_mass_flowrate_across_each_parallel_branch(
                    pressure_change_across_each_branch
                    );

            let mut mass_flowrate: MassRate = 
                mass_flowrate_across_each_branch[0];


            // get absolute value
            mass_flowrate = mass_flowrate.abs();

            mass_flowrate

        }
    // fluid mechanics calcs
    // now in a closure
    // I should probably code this as a function instead
    fn dracs_fluid_mechanics_calc_mass_rate(
        pipe_34: &InsulatedFluidComponent,
        pipe_33: &InsulatedFluidComponent,
        pipe_32: &InsulatedFluidComponent,
        pipe_31a: &InsulatedFluidComponent,
        static_mixer_61_label_31: &InsulatedFluidComponent,
        dhx_tube_side_30b: &NonInsulatedFluidComponent,
        dhx_tube_side_heat_exchanger_30: &NonInsulatedFluidComponent,
        dhx_tube_side_30a: &NonInsulatedFluidComponent,
        tchx_35a: &NonInsulatedFluidComponent,
        tchx_35b: &NonInsulatedFluidComponent,
        static_mixer_60_label_36: &InsulatedFluidComponent,
        pipe_36a: &InsulatedFluidComponent,
        pipe_37: &InsulatedFluidComponent,
        flowmeter_60_37a: &NonInsulatedFluidComponent,
        pipe_38: &InsulatedFluidComponent,
        pipe_39: &InsulatedFluidComponent,
    )-> MassRate {

        let mut dracs_hot_branch = 
            FluidComponentCollection::new_series_component_collection();

        dracs_hot_branch.clone_and_add_component(pipe_34);
        dracs_hot_branch.clone_and_add_component(pipe_33);
        dracs_hot_branch.clone_and_add_component(pipe_32);
        dracs_hot_branch.clone_and_add_component(pipe_31a);
        dracs_hot_branch.clone_and_add_component(static_mixer_61_label_31);
        dracs_hot_branch.clone_and_add_component(dhx_tube_side_30b);
        dracs_hot_branch.clone_and_add_component(dhx_tube_side_heat_exchanger_30);
        dracs_hot_branch.clone_and_add_component(dhx_tube_side_30a);


        let mut dracs_cold_branch = 
            FluidComponentCollection::new_series_component_collection();

        dracs_cold_branch.clone_and_add_component(tchx_35a);
        dracs_cold_branch.clone_and_add_component(tchx_35b);
        dracs_cold_branch.clone_and_add_component(static_mixer_60_label_36);
        dracs_cold_branch.clone_and_add_component(pipe_36a);
        dracs_cold_branch.clone_and_add_component(pipe_37);
        dracs_cold_branch.clone_and_add_component(flowmeter_60_37a);
        dracs_cold_branch.clone_and_add_component(pipe_38);
        dracs_cold_branch.clone_and_add_component(pipe_39);

        let mut dracs_branches = 
            FluidComponentSuperCollection::default();

        dracs_branches.set_orientation_to_parallel();
        dracs_branches.fluid_component_super_vector.push(dracs_hot_branch);
        dracs_branches.fluid_component_super_vector.push(dracs_cold_branch);

        let mass_rate = get_dracs_flowrate(&dracs_branches);

        mass_rate

    }

    // now the thermal hydraulics bit 
    fn calculate_dracs_thermal_hydraulics(
        mass_flowrate_counter_clockwise: MassRate,
        heat_rate_through_dhx: Power,
        tchx_heat_transfer_coeff: HeatTransfer,
        average_temperature_for_density_calcs: ThermodynamicTemperature,
        timestep: Time,
        pipe_34: &mut InsulatedFluidComponent,
        pipe_33: &mut InsulatedFluidComponent,
        pipe_32: &mut InsulatedFluidComponent,
        pipe_31a: &mut InsulatedFluidComponent,
        static_mixer_61_label_31: &mut InsulatedFluidComponent,
        dhx_tube_side_30b: &mut NonInsulatedFluidComponent,
        dhx_tube_side_heat_exchanger_30: &mut NonInsulatedFluidComponent,
        dhx_tube_side_30a: &mut NonInsulatedFluidComponent,
        tchx_35a: &mut NonInsulatedFluidComponent,
        tchx_35b: &mut NonInsulatedFluidComponent,
        static_mixer_60_label_36: &mut InsulatedFluidComponent,
        pipe_36a: &mut InsulatedFluidComponent,
        pipe_37: &mut InsulatedFluidComponent,
        flowmeter_60_37a: &mut NonInsulatedFluidComponent,
        pipe_38: &mut InsulatedFluidComponent,
        pipe_39: &mut InsulatedFluidComponent,
        ){

        // for an ideal situation, we have zero parasitic heat losses
        // therefore, for each component, except tchx, heat transfer 
        // coeff is zero

        let adiabatic_heat_transfer_coeff = HeatTransfer::ZERO;

        // create the heat transfer interaction 
        let advection_heat_transfer_interaction: HeatTransferInteractionType;

        // I'm going to create the advection interaction

        let average_therminol_density = 
            LiquidMaterial::TherminolVP1.density(
                average_temperature_for_density_calcs).unwrap();

        advection_heat_transfer_interaction = 
            HeatTransferInteractionType::
            new_advection_interaction(mass_flowrate_counter_clockwise, 
                                      average_therminol_density, 
                                      average_therminol_density);

        // now, let's link the fluid arrays using advection 
        // (no conduction here axially between arrays)
        {
            dhx_tube_side_30a.pipe_fluid_array.link_to_front(
                &mut dhx_tube_side_heat_exchanger_30.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            dhx_tube_side_heat_exchanger_30.pipe_fluid_array.link_to_front(
                &mut dhx_tube_side_30b.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            dhx_tube_side_30b.pipe_fluid_array.link_to_front(
                &mut static_mixer_61_label_31.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            static_mixer_61_label_31.pipe_fluid_array.link_to_front(
                &mut pipe_31a.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_31a.pipe_fluid_array.link_to_front(
                &mut pipe_32.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_32.pipe_fluid_array.link_to_front(
                &mut pipe_33.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_33.pipe_fluid_array.link_to_front(
                &mut pipe_34.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_34.pipe_fluid_array.link_to_front(
                &mut tchx_35a.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            tchx_35a.pipe_fluid_array.link_to_front(
                &mut tchx_35b.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            tchx_35b.pipe_fluid_array.link_to_front(
                &mut static_mixer_60_label_36.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            static_mixer_60_label_36.pipe_fluid_array.link_to_front(
                &mut pipe_36a.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_36a.pipe_fluid_array.link_to_front(
                &mut pipe_37.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_37.pipe_fluid_array.link_to_front(
                &mut flowmeter_60_37a.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            flowmeter_60_37a.pipe_fluid_array.link_to_front(
                &mut pipe_38.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_38.pipe_fluid_array.link_to_front(
                &mut pipe_39.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

            pipe_39.pipe_fluid_array.link_to_front(
                &mut dhx_tube_side_30a.pipe_fluid_array, 
                advection_heat_transfer_interaction)
                .unwrap();

        }
        // set the relevant heat transfer coefficients 
        // all zero except for tchx
        {
            // hot branch
            dhx_tube_side_30a.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            dhx_tube_side_heat_exchanger_30.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            dhx_tube_side_30b.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;

            static_mixer_61_label_31.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_31a.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;

            pipe_32.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_33.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_34.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;

            // cold branch 
            tchx_35a.heat_transfer_to_ambient = 
                tchx_heat_transfer_coeff;
            tchx_35b.heat_transfer_to_ambient = 
                tchx_heat_transfer_coeff;

            static_mixer_60_label_36.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_36a.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            flowmeter_60_37a.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;

            pipe_37.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_38.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;
            pipe_39.heat_transfer_to_ambient = 
                adiabatic_heat_transfer_coeff;

        }
        // add lateral heat losses and power through dhx
        {
            let zero_power: Power = Power::ZERO;

            // hot branch
            //
            // we add heat in through dhx 30 
            // everywhere else is zero heater power
            dhx_tube_side_30a
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            dhx_tube_side_heat_exchanger_30
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                heat_rate_through_dhx)
                .unwrap();
            dhx_tube_side_30b
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();

            static_mixer_61_label_31
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_31a
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_32
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_33
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_34
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();

            // cold branch 
            // ambient temperature of tchx is 20C  
            tchx_35a.ambient_temperature = 
                ThermodynamicTemperature::new::<degree_celsius>(20.0);
            tchx_35b.ambient_temperature = 
                ThermodynamicTemperature::new::<degree_celsius>(20.0);

            tchx_35a
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            tchx_35b
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();


            static_mixer_60_label_36
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_36a
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();

            flowmeter_60_37a.
                lateral_and_miscellaneous_connections_no_wall_correction(
                    mass_flowrate_counter_clockwise, 
                    zero_power)
                .unwrap();

            pipe_37
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_38
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();
            pipe_39
                .lateral_and_miscellaneous_connections_no_wall_correction(
                mass_flowrate_counter_clockwise, 
                zero_power)
                .unwrap();

        }
        
        // now we should be ready to advance timestep
        {
            dhx_tube_side_30a
                .advance_timestep(timestep)
                .unwrap();
            dhx_tube_side_heat_exchanger_30
                .advance_timestep(timestep)
                .unwrap();
            dhx_tube_side_30b
                .advance_timestep(timestep)
                .unwrap();

            static_mixer_61_label_31
                .advance_timestep(timestep)
                .unwrap();
            pipe_31a
                .advance_timestep(timestep)
                .unwrap();

            pipe_32
                .advance_timestep(timestep)
                .unwrap();
            pipe_33
                .advance_timestep(timestep)
                .unwrap();
            pipe_34
                .advance_timestep(timestep)
                .unwrap();

            // cold branch 
            tchx_35a
                .advance_timestep(timestep)
                .unwrap();
            tchx_35b
                .advance_timestep(timestep)
                .unwrap();

            static_mixer_60_label_36
                .advance_timestep(timestep)
                .unwrap();
            pipe_36a
                .advance_timestep(timestep)
                .unwrap();

            flowmeter_60_37a
                .advance_timestep(timestep)
                .unwrap();

            pipe_37
                .advance_timestep(timestep)
                .unwrap();
            pipe_38
                .advance_timestep(timestep)
                .unwrap();
            pipe_39
                .advance_timestep(timestep)
                .unwrap();
        }

        // we do it in serial, so it keeps things simple 
        // now we are done

    }
    
    // I also want to find the final temperature, which should be 
    // around the set point within thermocouple error (+/- 0.5 K)
    let mut final_tchx_outlet_temperature: ThermodynamicTemperature 
        = ThermodynamicTemperature::ZERO;

    // main simulation loop
    while current_simulation_time < max_simulation_time {
        // show the outlet temperature of tchx 

        let tchx_outlet_temperature: ThermodynamicTemperature = {

            // the front of the tchx is connected to static mixer 
            // 60 label 36
            let tchx35b_pipe_fluid_array_clone: FluidArray = 
                tchx_35b.pipe_fluid_array
                .clone()
                .try_into()
                .unwrap();

            // take the front single cv temperature 
            //
            // front single cv temperature is defunct
            // probably need to debug this

            let tchx_35b_front_single_cv_temperature: ThermodynamicTemperature 
                = tchx35b_pipe_fluid_array_clone
                .front_single_cv
                .temperature;



            let _tchx_35b_array_temperature: Vec<ThermodynamicTemperature>
                = tchx_35b
                .pipe_fluid_array_temperature()
                .unwrap();

            //dbg!(&tchx_35b_array_temperature);

            tchx_35b_front_single_cv_temperature

        };

        // we will need to change the tchx heat transfer coefficient 
        // using the PID controller
        //
        // record tchx outlet temperature if it is last 5s of time 
        
        let tchx_temperature_record_time_threshold = max_simulation_time - 
            Time::new::<second>(5.0);

        if current_simulation_time > tchx_temperature_record_time_threshold {
            final_tchx_outlet_temperature = tchx_outlet_temperature;
        }


        
        tchx_heat_transfer_coeff = {
            // first, calculate the set point error 

            let reference_temperature_interval_deg_celsius = 80.0;

            // error = y_sp - y_measured
            let set_point_abs_error_deg_celsius = 
                - tchx_outlet_temperature_set_point.get::<kelvin>()
                + tchx_outlet_temperature.get::<kelvin>();

            let nondimensional_error: Ratio = 
                (set_point_abs_error_deg_celsius/
                reference_temperature_interval_deg_celsius).into();

            // let's get the output 

            let dimensionless_heat_trf_input: Ratio
                = pid_controller.set_user_input_and_calc(
                nondimensional_error, 
                current_simulation_time).unwrap();

            // the dimensionless output is:
            //
            // (desired output - ref_val)/ref_val = dimensionless_input
            // 
            //
            // the reference value is decided by the user 
            // in this case 250 W/(m^2 K)
            
            let mut tchx_heat_trf_output = 
                dimensionless_heat_trf_input * reference_tchx_htc
                + reference_tchx_htc;

            // make sure it cannot be less than a certain amount 
            let tchx_minimum_heat_transfer = 
                HeatTransfer::new::<watt_per_square_meter_kelvin>(
                    10.0);

            // this makes it physically realistic
            if tchx_heat_trf_output < tchx_minimum_heat_transfer {
                tchx_heat_trf_output = tchx_minimum_heat_transfer;
            }

            tchx_heat_trf_output

        };
        // fluid first 
        // 
        let mass_flowrate_absolute: MassRate = 
            dracs_fluid_mechanics_calc_mass_rate(
                &pipe_34, 
                &pipe_33, 
                &pipe_32, 
                &pipe_31a, 
                &static_mixer_61_label_31, 
                &dhx_tube_side_30b, 
                &dhx_tube_side_heat_exchanger_30, 
                &dhx_tube_side_30a, 
                &tchx_35a, 
                &tchx_35b, 
                &static_mixer_60_label_36, 
                &pipe_36a, 
                &pipe_37, 
                &flowmeter_60_37a, 
                &pipe_38, 
                &pipe_39);


        // I assume the mass_flowrate_counter_clockwise 
        // is the same as absolute flowrate
        let mass_flowrate_counter_clockwise = 
            mass_flowrate_absolute;

        // next, thermal hydraulics calcs 

        calculate_dracs_thermal_hydraulics(
            mass_flowrate_counter_clockwise, 
            heat_rate_through_dhx, 
            tchx_heat_transfer_coeff, 
            average_temperature_for_density_calcs, 
            timestep, 
            &mut pipe_34, 
            &mut pipe_33, 
            &mut pipe_32, 
            &mut pipe_31a, 
            &mut static_mixer_61_label_31, 
            &mut dhx_tube_side_30b, 
            &mut dhx_tube_side_heat_exchanger_30, 
            &mut dhx_tube_side_30a, 
            &mut tchx_35a, 
            &mut tchx_35b, 
            &mut static_mixer_60_label_36, 
            &mut pipe_36a, 
            &mut pipe_37, 
            &mut flowmeter_60_37a, 
            &mut pipe_38, 
            &mut pipe_39);






        current_simulation_time += timestep;
        let debug: bool = false;
        if debug {
            // show the mass flowrate
            // tchx outlet temperature 
            // current sim time 
            // and tchx heat trf coeff
            dbg!(&mass_flowrate_absolute);
            dbg!(&tchx_outlet_temperature);
            dbg!(&current_simulation_time);
            dbg!(&tchx_heat_transfer_coeff);
        }

    }

    // panic to see debug messages

    //panic!();
    
    // final assertion 
    // that tchx outlet temperature is equal to set point within 
    // thermocouple error 

    let thermocouple_error_kelvin: f64 = 0.5;

    approx::assert_abs_diff_eq!(
        tchx_outlet_temperature_set_point.get::<degree_celsius>(),
        final_tchx_outlet_temperature.get::<degree_celsius>(),
        epsilon=thermocouple_error_kelvin
        );


}


