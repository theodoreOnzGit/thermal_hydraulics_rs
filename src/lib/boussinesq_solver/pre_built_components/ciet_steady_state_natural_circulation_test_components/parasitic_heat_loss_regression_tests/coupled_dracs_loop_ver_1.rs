#[test]
pub fn test_uncalibrated_dracs_loop(){

    verify_coupled_dhx_analytical_solution_version_1(
        2764.53, 
        40.0,
        4.6990e-2,
        3.4570e-2,
        ).unwrap();
}

#[cfg(test)]
/// function to verify the dhx analytical solution
/// version 1,
/// the DHX here uses uncalibrated Gnielinski correlations 
/// to estimate heat transfer coefficients
pub fn verify_coupled_dhx_analytical_solution_version_1(
    input_power_watts: f64,
    tchx_outlet_temperature_set_point_degc: f64,
    experimental_dracs_mass_flowrate_kg_per_s: f64,
    experimental_primary_mass_flowrate_kg_per_s: f64) -> 
Result<(),crate::thermal_hydraulics_error::ThermalHydraulicsLibError>{
    use uom::si::{f64::*, mass_rate::kilogram_per_second, power::watt};

    use uom::si::{frequency::hertz, ratio::ratio, time::millisecond};

    use crate::boussinesq_solver::pre_built_components::ciet_isothermal_test_components::*;
    use crate::boussinesq_solver::pre_built_components::ciet_steady_state_natural_circulation_test_components::coupled_dracs_loop_tests::dhx_constructor::new_dhx_sthe_version_1;
    use crate::boussinesq_solver::pre_built_components::ciet_steady_state_natural_circulation_test_components::coupled_dracs_loop_tests::dracs_loop_calc_functions::{coupled_dracs_fluid_mechanics_calc_abs_mass_rate, coupled_dracs_loop_link_up_components, dracs_loop_advance_timestep_except_dhx};
    use crate::boussinesq_solver::pre_built_components::ciet_steady_state_natural_circulation_test_components::coupled_dracs_loop_tests::pri_loop_calc_functions::{coupled_dracs_pri_loop_branches_fluid_mechanics_calc_abs_mass_rate, coupled_dracs_pri_loop_dhx_heater_link_up_components, pri_loop_advance_timestep_except_dhx, pri_loop_dhx_shell_temperature_diagnostics, pri_loop_heater_temperature_diagnostics};
    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::dracs_loop_components::*;
    use crate::prelude::beta_testing::FluidArray;
    use uom::ConstZero;

    use uom::si::thermodynamic_temperature::{degree_celsius, kelvin};
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use uom::si::time::second;

    let input_power = Power::new::<watt>(input_power_watts);
    let experimental_dracs_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_dracs_mass_flowrate_kg_per_s);
    let experimental_primary_mass_flowrate = 
        MassRate::new::<kilogram_per_second>(
            experimental_primary_mass_flowrate_kg_per_s);

    let tchx_outlet_temperature_set_point = 
        ThermodynamicTemperature::new::<degree_celsius>(
            tchx_outlet_temperature_set_point_degc);
    use chem_eng_real_time_process_control_simulator::alpha_nightly::transfer_fn_wrapper_and_enums::TransferFnTraits;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::ProportionalController;
    use chem_eng_real_time_process_control_simulator::alpha_nightly::controllers::AnalogController;

    // max error is 0.5% according to SAM 
    // is okay, because typical flowmeter measurement error is 2% anyway
    let timestep = Time::new::<second>(0.05);
    let heat_rate_through_heater = input_power;
    let mut tchx_heat_transfer_coeff: HeatTransfer;

    let reference_tchx_htc = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(40.0);
    let average_temperature_for_density_calcs = 
        ThermodynamicTemperature::new::<degree_celsius>(80.0);
    // let's calculate 400 seconds of simulated time 
    // it takes about that long for the temperature to settle down
    // this is compared to value at 4000s

    let mut current_simulation_time = Time::ZERO;
    let max_simulation_time = Time::new::<second>(400.0);

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



    let initial_temperature = tchx_outlet_temperature_set_point;

    // DRACS hot branch or (mostly) hot leg
    let mut pipe_34 = new_pipe_34(initial_temperature);
    let mut pipe_33 = new_pipe_33(initial_temperature);
    let mut pipe_32 = new_pipe_32(initial_temperature);
    let mut pipe_31a = new_pipe_31a(initial_temperature);
    let mut static_mixer_61_label_31 = new_static_mixer_61_label_31(initial_temperature);
    let mut dhx_tube_side_30b = new_dhx_tube_side_30b(initial_temperature);
    let mut dhx_sthe = new_dhx_sthe_version_1(initial_temperature);
    let mut dhx_tube_side_30a = new_dhx_tube_side_30a(initial_temperature);

    // DRACS cold branch or (mostly) cold leg
    let mut tchx_35a = new_ndhx_tchx_horizontal_35a(initial_temperature);
    let mut tchx_35b = new_ndhx_tchx_vertical_35b(initial_temperature);
    let mut static_mixer_60_label_36 = new_static_mixer_60_label_36(initial_temperature);
    let mut pipe_36a = new_pipe_36a(initial_temperature);
    let mut pipe_37 = new_pipe_37(initial_temperature);
    let mut flowmeter_60_37a = new_flowmeter_60_37a(initial_temperature);
    let mut pipe_38 = new_pipe_38(initial_temperature);
    let mut pipe_39 = new_pipe_39(initial_temperature);

    // pri loop dhx branch top to bottom 5a to 17b 

    let mut pipe_5a = new_branch_5a(initial_temperature);
    let mut pipe_26 = new_pipe_26(initial_temperature);
    let mut pipe_25a = new_pipe_25a(initial_temperature);
    let mut static_mixer_21_label_25 = new_static_mixer_21_label_25(initial_temperature);
    // here is where the dhx shell side should be (component 24)
    let mut pipe_23a = new_pipe_23a(initial_temperature);
    let mut static_mixer_20_label_23 = new_static_mixer_20_label_23(initial_temperature);
    let mut pipe_22 = new_pipe_22(initial_temperature);
    let mut flowmeter_20_21a = new_flowmeter_20_label_21a(initial_temperature);
    let mut pipe_21 = new_pipe_21(initial_temperature);
    let mut pipe_20 = new_pipe_20(initial_temperature);
    let mut pipe_19 = new_pipe_19(initial_temperature);
    let mut pipe_17b = new_branch_17b(initial_temperature);

    // heater branch top to bottom 4 to 18
    let mut pipe_4 = new_pipe_4(initial_temperature);
    let mut pipe_3 = new_pipe_3(initial_temperature);
    let mut pipe_2a = new_pipe_2a(initial_temperature);
    let mut static_mixer_10_label_2 = new_static_mixer_10_label_2(initial_temperature);
    let mut heater_top_head_1a = new_heater_top_head_1a(initial_temperature);
    let mut heater_ver_1 = new_heated_section_version_1_label_1(initial_temperature);
    let mut heater_bottom_head_1b = new_heater_bottom_head_1b(initial_temperature);
    let mut pipe_18 = new_pipe_18(initial_temperature);


    let mut final_mass_flowrate_pri_loop: MassRate 
        = MassRate::ZERO;
    let mut final_mass_flowrate_dracs_loop: MassRate 
        = MassRate::ZERO;
    let mut _final_tchx_outlet_temperature: ThermodynamicTemperature 
        = ThermodynamicTemperature::ZERO;

    let ambient_htc = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // calculation loop
    while current_simulation_time < max_simulation_time {

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
            _final_tchx_outlet_temperature = tchx_outlet_temperature;
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
                    5.0);

            // this makes it physically realistic
            if tchx_heat_trf_output < tchx_minimum_heat_transfer {
                tchx_heat_trf_output = tchx_minimum_heat_transfer;
            }

            tchx_heat_trf_output

        };

        // fluid calculation loop 
        //
        // first, absolute mass flowrate across two branches
        let dhx_tube_side_heat_exchanger_30 = 
            dhx_sthe.get_clone_of_tube_side_parallel_tube_fluid_component();
        let dhx_shell_side_pipe_24 = 
            dhx_sthe.get_clone_of_shell_side_fluid_component();

        let absolute_mass_flowrate_dracs = 
            coupled_dracs_fluid_mechanics_calc_abs_mass_rate(
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

        // likely the natural circulation is counter clockwise 
        let counter_clockwise_dracs_flowrate = absolute_mass_flowrate_dracs;

        let absolute_mass_flowrate_pri_loop = 
            coupled_dracs_pri_loop_branches_fluid_mechanics_calc_abs_mass_rate(
                &pipe_4, 
                &pipe_3, 
                &pipe_2a, 
                &static_mixer_10_label_2, 
                &heater_top_head_1a, 
                &heater_ver_1, 
                &heater_bottom_head_1b, 
                &pipe_18, 
                &pipe_5a, 
                &pipe_26, 
                &pipe_25a, 
                &static_mixer_21_label_25, 
                &dhx_shell_side_pipe_24, 
                &static_mixer_20_label_23, 
                &pipe_23a, 
                &pipe_22, 
                &flowmeter_20_21a, 
                &pipe_21, 
                &pipe_20, 
                &pipe_19, 
                &pipe_17b);

        let counter_clockwise_pri_loop_flowrate = absolute_mass_flowrate_pri_loop;

        // next, link up the heat transfer entities 
        // all lateral linking is done except for DHX
        coupled_dracs_loop_link_up_components(
            counter_clockwise_dracs_flowrate, 
            tchx_heat_transfer_coeff, 
            average_temperature_for_density_calcs, 
            timestep, 
            ambient_htc, 
            &mut pipe_34, 
            &mut pipe_33, 
            &mut pipe_32, 
            &mut pipe_31a, 
            &mut static_mixer_61_label_31, 
            &mut dhx_tube_side_30b, 
            &mut dhx_sthe, 
            &mut dhx_tube_side_30a, 
            &mut tchx_35a, 
            &mut tchx_35b, 
            &mut static_mixer_60_label_36, 
            &mut pipe_36a, 
            &mut pipe_37, 
            &mut flowmeter_60_37a, 
            &mut pipe_38, 
            &mut pipe_39);

        coupled_dracs_pri_loop_dhx_heater_link_up_components(
            counter_clockwise_pri_loop_flowrate, 
            heat_rate_through_heater, 
            average_temperature_for_density_calcs, 
            ambient_htc, 
            &mut pipe_4, 
            &mut pipe_3, 
            &mut pipe_2a, 
            &mut static_mixer_10_label_2, 
            &mut heater_top_head_1a, 
            &mut heater_ver_1, 
            &mut heater_bottom_head_1b, 
            &mut pipe_18, 
            &mut pipe_5a, 
            &mut pipe_26, 
            &mut pipe_25a, 
            &mut static_mixer_21_label_25, 
            &mut dhx_sthe, 
            &mut static_mixer_20_label_23, 
            &mut pipe_23a, 
            &mut pipe_22, 
            &mut flowmeter_20_21a, 
            &mut pipe_21, 
            &mut pipe_20, 
            &mut pipe_19, 
            &mut pipe_17b);

        // advance timestep
        dracs_loop_advance_timestep_except_dhx(
            timestep, &mut pipe_34, &mut pipe_33, &mut pipe_32, 
            &mut pipe_31a, &mut static_mixer_61_label_31, 
            &mut dhx_tube_side_30b, &mut dhx_tube_side_30a, 
            &mut tchx_35a, &mut tchx_35b, &mut static_mixer_60_label_36, 
            &mut pipe_36a, &mut pipe_37, &mut flowmeter_60_37a, 
            &mut pipe_38, &mut pipe_39);

        pri_loop_advance_timestep_except_dhx(
            timestep, &mut pipe_4, &mut pipe_3, &mut pipe_2a, 
            &mut static_mixer_10_label_2, &mut heater_top_head_1a, 
            &mut heater_ver_1, &mut heater_bottom_head_1b, 
            &mut pipe_18, &mut pipe_5a, &mut pipe_26, &mut pipe_25a, 
            &mut static_mixer_21_label_25, &mut static_mixer_20_label_23, 
            &mut pipe_23a, &mut pipe_22, &mut flowmeter_20_21a, 
            &mut pipe_21, &mut pipe_20, &mut pipe_19, &mut pipe_17b);

        // for dhx, a little more care is needed to do the 
        // lateral and misc connections and advance timestep 
        // advance timestep
        //
        // by default, dhx flowrate is downwards in this setup

        let prandtl_wall_correction_setting = true; 
        let tube_side_total_mass_flowrate = -counter_clockwise_dracs_flowrate;
        let shell_side_total_mass_flowrate = counter_clockwise_pri_loop_flowrate;

        dhx_sthe.lateral_and_miscellaneous_connections(
            prandtl_wall_correction_setting, 
            tube_side_total_mass_flowrate, 
            shell_side_total_mass_flowrate).unwrap();

        dhx_sthe.advance_timestep(timestep).unwrap();

        

        // record 
        if current_simulation_time > tchx_temperature_record_time_threshold {
            final_mass_flowrate_dracs_loop = counter_clockwise_dracs_flowrate;
            final_mass_flowrate_pri_loop = counter_clockwise_pri_loop_flowrate;
        }

        // debugging 
        let debug_settings = true;

        // temperatures before and after heater
        let ((_bt_11,_wt_10),(_bt_12,_wt_13)) = 
            pri_loop_heater_temperature_diagnostics(
            &mut heater_bottom_head_1b, 
            &mut static_mixer_10_label_2, 
            debug_settings);
        // temperatures before and after dhx shell
        let ((_bt_21,_wt_20),(_bt_27,_wt_26)) = 
            pri_loop_dhx_shell_temperature_diagnostics(
            &mut pipe_25a, 
            &mut static_mixer_20_label_23, 
            debug_settings);
        

        current_simulation_time += timestep;

    }

    dbg!(&(
            input_power,
            final_mass_flowrate_pri_loop,
            final_mass_flowrate_dracs_loop
            ));

    approx::assert_relative_eq!(
        experimental_primary_mass_flowrate.get::<kilogram_per_second>(),
        final_mass_flowrate_pri_loop.get::<kilogram_per_second>(),
        max_relative=0.01);

    approx::assert_relative_eq!(
        experimental_dracs_mass_flowrate.get::<kilogram_per_second>(),
        final_mass_flowrate_dracs_loop.get::<kilogram_per_second>(),
        max_relative=0.01);

    Ok(())

}
