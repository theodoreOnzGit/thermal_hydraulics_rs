


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
    use crate::boussinesq_solver::pre_built_components::
        ciet_steady_state_natural_circulation_test_components::
        dracs_loop_components::new_isolated_dhx_tube_side_30;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use std::time::SystemTime;

    use uom::si::power::kilowatt;
    use uom::ConstZero;
    use uom::si::time::second;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::heat_transfer::watt_per_square_meter_kelvin;

    use crate::boussinesq_solver::boundary_conditions::BCType;
    use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
    use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
    use crate::boussinesq_solver::pre_built_components::heat_transfer_entities::HeatTransferEntity;

    let initial_temperature = 
        ThermodynamicTemperature::new::<degree_celsius>(25.0);

    let mut adiabatic_dhx_tube_side_30 = 
        new_isolated_dhx_tube_side_30(initial_temperature);

    let htc_to_ambient_zero = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(0.0);

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

    // main loop
    
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
    dbg!(&adiabatic_dhx_tube_side_30_outlet_temp);
    todo!()

}
