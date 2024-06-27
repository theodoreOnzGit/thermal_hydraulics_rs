





// for coupled dracs loop, first thing first is to 
// debug parallel bare pipes 
//
// first thing first, heat addition 
//
// suppose we have 9 kW of heat into one pipe 
// at 0.18 kg/s TherminolVP1
//
// whether you have the flow split into 1 single pipe or 10 parallel 
// pipes should not matter. After a set time, the temperature should 
// be the same. presuming of course, no heat loss to environment
//
// I'm taking dhx tube side 30b as a reference
// but setting heat transfer to ambient as 0
//
#[test]
pub fn parallel_bare_pipes_debugging_heat_addition(){

    use uom::si::f64::*;
    use uom::si::ratio::ratio;
    use uom::si::length::{meter, millimeter};
    use uom::si::area::square_meter;
    use uom::si::angle::degree;

    use crate::boussinesq_solver::{array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray, boussinesq_thermophysical_properties::SolidMaterial, heat_transfer_correlations::nusselt_number_correlations::enums::NusseltCorrelation, pre_built_components::non_insulated_fluid_components::NonInsulatedFluidComponent};
    use uom::si::thermodynamic_temperature::degree_celsius;

    use uom::si::power::kilowatt;
    use uom::ConstZero;
    use uom::si::time::second;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;
    use uom::si::pressure::atmosphere;

    use crate::boussinesq_solver::boundary_conditions::BCType;
    use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
    use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
    use crate::boussinesq_solver::pre_built_components::heat_transfer_entities::HeatTransferEntity;

    use crate::boussinesq_solver::pre_built_components::
        non_insulated_parallel_fluid_components::
        NonInsulatedParallelFluidComponent;
    // conditions 
    let initial_temperature = 
        ThermodynamicTemperature::new::<degree_celsius>(25.0);


    let htc_to_ambient_zero = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(0.0);


    // parameters for new_isolated_dhx_tube_side_30
    let ambient_temperature = ThermodynamicTemperature::new::<degree_celsius>(20.0);
    let fluid_pressure = Pressure::new::<atmosphere>(1.0);
    let solid_pressure = Pressure::new::<atmosphere>(1.0);
    // for dhx modelling,
    //
    // dhx tubes are modelled in SAM as 19 tubes of diameter 
    // 0.00635 m 
    // and flow area of 6.1072e-4 m^2
    //
    // in Zweibaum's RELAP model,
    // it is quite different from the SAM model 
    // which follows Bickel's data
    let hydraulic_diameter = Length::new::<meter>(6.35e-3);
    let pipe_length = Length::new::<meter>(1.18745);
    let flow_area = Area::new::<square_meter>(6.0172e-4);
    let incline_angle = Angle::new::<degree>(90.0-180.0);
    let form_loss = Ratio::new::<ratio>(3.3);
    //estimated component wall roughness (doesn't matter here,
    //but i need to fill in)
    let surface_roughness = Length::new::<millimeter>(0.015);
    let id = hydraulic_diameter;
    let pipe_thickness = Length::new::<meter>(0.00079375);
    let od = id + pipe_thickness;
    let pipe_shell_material = SolidMaterial::SteelSS304L;
    let pipe_fluid = LiquidMaterial::TherminolVP1;
    let htc_to_ambient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // from SAM nodalisation, we have 11 nodes only, 
    // now because there are two outer nodes, the 
    // number of inner nodes is 11-2
    let user_specified_inner_nodes = 11-2; 



    let mut adiabatic_dhx_tube_side_30 = 
        NonInsulatedFluidComponent::new_bare_pipe(
            initial_temperature, 
            ambient_temperature, 
            fluid_pressure, 
            solid_pressure, 
            flow_area, 
            incline_angle, 
            form_loss, 
            id, 
            od, 
            pipe_length, 
            hydraulic_diameter, 
            surface_roughness, 
            pipe_shell_material, 
            pipe_fluid, 
            htc_to_ambient, 
            user_specified_inner_nodes);

    // for heat exchangers, I give an ideal Nusselt number correlation 
    // as an approximation so that film thermal resistance is minimised
    let mut fluid_array_ideal_nusslet: FluidArray = 
        adiabatic_dhx_tube_side_30.pipe_fluid_array
        .clone()
        .try_into()
        .unwrap();

    fluid_array_ideal_nusslet.nusselt_correlation = 
        NusseltCorrelation::IdealNusseltOneBillion;

    adiabatic_dhx_tube_side_30.pipe_fluid_array = 
        fluid_array_ideal_nusslet.into();

    adiabatic_dhx_tube_side_30.heat_transfer_to_ambient = 
        htc_to_ambient_zero;

    // now let's do a simple loop to check temperature after 30s
    let max_time = Time::new::<second>(10.0);
    let timestep = Time::new::<second>(0.01);
    let mut simulation_time = Time::ZERO;
    let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);
    let heater_power = Power::new::<kilowatt>(9.0);

    let mut inlet_bc: HeatTransferEntity = BCType::new_const_temperature( 
        initial_temperature).into();

    let mut outlet_bc: HeatTransferEntity = BCType::new_adiabatic_bc().into();

    let mut adiabatic_dhx_tube_side_30_outlet_temp = 
        ThermodynamicTemperature::ZERO;

    // main loop for single tube
    
    while max_time > simulation_time {

        // get average temperature and advection first 
        // we use boussinesq approximation, so the average temperature 
        // for density doesn't change
        // i'll just arbitrarily do 50 C
        let average_temperature_for_density_calcs = 
            ThermodynamicTemperature::
            new::<degree_celsius>(50.0);

        let average_therminol_density = 
            LiquidMaterial::TherminolVP1.density(
                average_temperature_for_density_calcs).unwrap();

        let advection_heat_transfer_interaction = 
            HeatTransferInteractionType::
            new_advection_interaction(mass_flowrate, 
                                      average_therminol_density, 
                                      average_therminol_density);

        // now let's link the fluid array part first 

        {
            adiabatic_dhx_tube_side_30.pipe_fluid_array.link_to_back(
                &mut inlet_bc, 
                advection_heat_transfer_interaction
                ).unwrap();

            adiabatic_dhx_tube_side_30.pipe_fluid_array.link_to_front(
                &mut outlet_bc, 
                advection_heat_transfer_interaction
                ).unwrap();
        }

        // set htc to zero (again, just kiasu)
        {
            adiabatic_dhx_tube_side_30.heat_transfer_to_ambient 
                = htc_to_ambient_zero;
        }
        // add 9 kW as a heater power for reference 
        {
            adiabatic_dhx_tube_side_30
                .lateral_and_miscellaneous_connections(
                    mass_flowrate, 
                    heater_power).unwrap();
        }

        // now advance timestep 
        {
            adiabatic_dhx_tube_side_30.
                advance_timestep(timestep).unwrap();
        }

        // let's obtain temperature 

        // first i get the array, then take the last element of the array
        // this is the "front_cv"
        let outlet_temperature: ThermodynamicTemperature 
            = *adiabatic_dhx_tube_side_30
            .pipe_fluid_array_temperature()
            .unwrap()
            .clone()
            .last()
            .unwrap();


        adiabatic_dhx_tube_side_30_outlet_temp = 
            outlet_temperature;



        
        // moving on timestep
        simulation_time += timestep;

    }
    // reset simulation time
    simulation_time = Time::ZERO;

    // suppose there is a parallel tube bundle of 20 equivalent 
    // dhx tube side 30
    let number_of_tubes = 20;
    let mass_flowrate_through_tube_bundle = 
        number_of_tubes as f64 * mass_flowrate;

    let heater_power_for_tube_bundle = 
        number_of_tubes as f64 * heater_power;

    let mut parallel_adiabatic_dhx_tube_side_30_outlet_temp = 
        ThermodynamicTemperature::ZERO;

    let mut parallel_adiabatic_dhx_tube_side_30 = 
        NonInsulatedParallelFluidComponent::new_bare_pipe_parallel_array(
            initial_temperature, 
            ambient_temperature, 
            fluid_pressure, 
            solid_pressure, 
            flow_area, 
            incline_angle, 
            form_loss, 
            id, 
            od, 
            pipe_length, 
            hydraulic_diameter, 
            surface_roughness, 
            pipe_shell_material, 
            pipe_fluid, 
            htc_to_ambient, 
            user_specified_inner_nodes,
            number_of_tubes);

    let mut fluid_array_ideal_nusslet: FluidArray = 
        adiabatic_dhx_tube_side_30.pipe_fluid_array
        .clone()
        .try_into()
        .unwrap();

    fluid_array_ideal_nusslet.nusselt_correlation = 
        NusseltCorrelation::IdealNusseltOneBillion;

    parallel_adiabatic_dhx_tube_side_30.pipe_fluid_array = 
        fluid_array_ideal_nusslet.into();

    parallel_adiabatic_dhx_tube_side_30.heat_transfer_to_ambient = 
        htc_to_ambient_zero;


    // main loop for parallel tube
    
    while max_time > simulation_time {

        // get average temperature and advection first 
        // we use boussinesq approximation, so the average temperature 
        // for density doesn't change
        // i'll just arbitrarily do 50 C
        let average_temperature_for_density_calcs = 
            ThermodynamicTemperature::
            new::<degree_celsius>(50.0);

        let average_therminol_density = 
            LiquidMaterial::TherminolVP1.density(
                average_temperature_for_density_calcs).unwrap();

        let advection_heat_transfer_interaction = 
            HeatTransferInteractionType::
            new_advection_interaction(mass_flowrate_through_tube_bundle, 
                                      average_therminol_density, 
                                      average_therminol_density);

        // now let's link the fluid array part first 

        {
            parallel_adiabatic_dhx_tube_side_30.pipe_fluid_array.link_to_back(
                &mut inlet_bc, 
                advection_heat_transfer_interaction
                ).unwrap();

            parallel_adiabatic_dhx_tube_side_30.pipe_fluid_array.link_to_front(
                &mut outlet_bc, 
                advection_heat_transfer_interaction
                ).unwrap();
        }

        // set htc to zero (again, just kiasu)
        {
            parallel_adiabatic_dhx_tube_side_30.heat_transfer_to_ambient 
                = htc_to_ambient_zero;
        }
        // add 9 kW as a heater power for reference 
        // times 20 for the bundle
        {
            parallel_adiabatic_dhx_tube_side_30
                .lateral_and_miscellaneous_connections(
                    mass_flowrate_through_tube_bundle, 
                    heater_power_for_tube_bundle).unwrap();
        }

        // now advance timestep 
        {
            parallel_adiabatic_dhx_tube_side_30.
                advance_timestep(timestep).unwrap();
        }

        // let's obtain temperature 

        // first i get the array, then take the last element of the array
        // this is the "front_cv"
        let outlet_temperature: ThermodynamicTemperature 
            = *parallel_adiabatic_dhx_tube_side_30
            .pipe_fluid_array_temperature()
            .unwrap()
            .clone()
            .last()
            .unwrap();


        parallel_adiabatic_dhx_tube_side_30_outlet_temp = 
            outlet_temperature;



        
        // moving on timestep
        simulation_time += timestep;

    }

    dbg!(&adiabatic_dhx_tube_side_30_outlet_temp);
    dbg!(&parallel_adiabatic_dhx_tube_side_30_outlet_temp);
    todo!()

}
