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
    SurfaceArea, SingleCVNode, CVType, BCType};
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
#[ignore = "takes about 20min, only use for data collection"]
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

    
    

    let therminol_mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);

    // construct the objects

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

    let inlet_const_temp_ptr = Arc::new(Mutex::new(
        inlet_const_temp
    ));

    let outlet_zero_heat_flux_ptr = Arc::new(Mutex::new(
        outlet_zero_heat_flux
    ));

    // the two types of HeatTransferInteractionType are 
    // advection and convection resistance
    //
    // 2007 square_centimeter
    // and 607 watt_per_square_meter_kelvin

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



    // timestep settings


    let max_time: Time = Time::new::<second>(20000.0);
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

            // todo calculation steps

            let mut therminol_cylinder_in_loop = 
            therminol_cylinder_ptr.lock().unwrap();
            let mut steel_shell_in_loop = 
            steel_shell_ptr.lock().unwrap();
            let mut inlet_const_temp_in_loop = 
            inlet_const_temp_ptr.lock().unwrap();
            let mut outlet_zero_heat_flux_in_loop = 
            outlet_zero_heat_flux_ptr.lock().unwrap();



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
