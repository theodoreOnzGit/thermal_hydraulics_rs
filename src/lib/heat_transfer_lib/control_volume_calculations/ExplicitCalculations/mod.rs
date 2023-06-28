extern crate uom;
use ndarray::FixedInitializer;
use uom::si::f64::*;
use uom::si::mass_rate::kilogram_per_second;
use crate::heat_transfer_lib::control_volume_calculations::Sandbox::*;
use uom::si::time::second;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::power::watt;

#[cfg(test)]
mod explicit_calc_sandbox {
    #[test]
    pub fn test_FluidEntityCollectionV1(){

        extern crate approx;

        use crate::heat_transfer_lib::ControlVolumeCalculations::Sandbox::
            v2_IterativeHeatFluxTherminolPipe;

        use crate::heat_transfer_lib::ControlVolumeCalculations::FluidEntity_StructsAndTraits::
            FluidEntityInitialisationSteps;

        use uom::si::f64::*;
        use uom::si::time::second;
        use uom::si::thermodynamic_temperature::kelvin;
        use uom::si::volume::cubic_meter;


        let timestep = Time::new::<second>(0.1_f64);
        let initial_global_temp = ThermodynamicTemperature::
            new::<kelvin>(300_f64);
        let fluid_volume = Volume::new::<cubic_meter>(
            0.01_f64.powf(3_f64));

        use crate::heat_transfer_lib::ControlVolumeCalculations::ExplicitCalculations::
            FluidEntityCollectionV1;

        // instantiate a new fluid entity collection object
        let mut fluid_entity_collection_obj = 
            FluidEntityCollectionV1::new();

        // initiate timestep and global initial temp

        fluid_entity_collection_obj.setup_step_0_set_timestep_and_initial_temp(
            timestep, 
            initial_global_temp);

        // then we add pipe 1, pipe 2 and pipe 3
        // First let's check the enthalpy at 300K
        //



        fluid_entity_collection_obj.setup_step_1_add_new_component(
            "pipe1".to_string(), 
            fluid_volume);

        // index 1
        fluid_entity_collection_obj.setup_step_1_add_new_component(
            "pipe2".to_string(), 
            fluid_volume);

        // index 2
        fluid_entity_collection_obj.setup_step_1_add_new_component(
            "pipe3".to_string(), 
            fluid_volume);


        // now to connect the pipes in this fashion
        // 1 -> 2 -> 3 
        // and 3 connects back to 1 in a circular fashion

        fluid_entity_collection_obj.
            setup_step_2_connect_inlet_and_outlet_pipe(0, 1);

        fluid_entity_collection_obj.
            setup_step_2_connect_inlet_and_outlet_pipe(1, 2);

        fluid_entity_collection_obj.
            setup_step_2_connect_inlet_and_outlet_pipe(2, 0);


        // if connected correctly, this should pass:
        let mut pipe1 = fluid_entity_collection_obj.
            fluid_entity_vector[0].clone();

        assert_eq!(0, pipe1.fluid_parameters.index_data.fluid_entity_index);

        assert_eq!(2, pipe1.fluid_parameters.index_data.
                   inlet_fluid_entity_index);

        assert_eq!(1, pipe1.fluid_parameters.index_data.
                   outlet_fluid_entity_index);

        let mut pipe2 = fluid_entity_collection_obj.
            fluid_entity_vector[1].clone();

        assert_eq!(1, pipe2.fluid_parameters.index_data.
                   fluid_entity_index);

        assert_eq!(0, pipe2.fluid_parameters.index_data.
                   inlet_fluid_entity_index);

        assert_eq!(2, pipe2.fluid_parameters.index_data.
                   outlet_fluid_entity_index);

        let mut pipe3 = fluid_entity_collection_obj.
            fluid_entity_vector[2].clone();

        assert_eq!(2, pipe3.fluid_parameters.index_data.
                   fluid_entity_index);
        assert_eq!(1, pipe3.fluid_parameters.index_data.
                   inlet_fluid_entity_index);
        assert_eq!(0, pipe3.fluid_parameters.index_data.
                   outlet_fluid_entity_index);


        // after adding and connecting
        // these pipes, i must have their enthalpy equal to
        // that at 300K

        use crate::ControlVolumeCalculations::TherminolDowthermPipes::
            TherminolFluidProperties;

        pub struct TestDowthermAProps {}
        impl TherminolFluidProperties for TestDowthermAProps {}

        let reference_enthalpy = TestDowthermAProps::
            enthalpy(initial_global_temp);

        // panic!("{}",reference_enthalpy.value.to_string());
        //
        let temp_300k = TestDowthermAProps::
            get_temperature_from_enthalpy(reference_enthalpy);

        approx::assert_relative_eq!(300.0,temp_300k.value,
                                    max_relative = 0.01);

        fluid_entity_collection_obj.
            step_1_calculate_current_timestep_temp_enthalpies();


        pipe1 = fluid_entity_collection_obj.
            fluid_entity_vector[0].clone();

        approx::assert_relative_eq!(
            300.0,
            pipe1.fluid_parameters.temperature_data.
            inlet_temp_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            300.0,
            pipe1.fluid_parameters.temperature_data.
            outlet_temp_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            300.0,
            pipe1.fluid_parameters.temperature_data.
            fluid_temp_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe1.fluid_parameters.enthalpy_data.
            inlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe1.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe1.fluid_parameters.enthalpy_data.
            fluid_enthalpy_old.value,
            max_relative = 0.001);


        pipe2 = fluid_entity_collection_obj.
            fluid_entity_vector[1].clone();

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe2.fluid_parameters.enthalpy_data.
            inlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe2.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe2.fluid_parameters.enthalpy_data.
            fluid_enthalpy_old.value,
            max_relative = 0.001);

        pipe3 = fluid_entity_collection_obj.
            fluid_entity_vector[2].clone();

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe3.fluid_parameters.enthalpy_data.
            inlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe3.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old.value,
            max_relative = 0.001);

        approx::assert_relative_eq!(
            reference_enthalpy.value,
            pipe3.fluid_parameters.enthalpy_data.
            fluid_enthalpy_old.value,
            max_relative = 0.001);

        // after adding these temperatures, we should assert that
        // the enthalpy is the same as the enthalpy at 300K 
        // which is the global temperature
        // next step is to set the mass flowrate, 
        // heat inputs and work inputs
        //

        use uom::si::mass_rate::kilogram_per_second;

        let mass_flowrate = 
            MassRate::new::<kilogram_per_second>(0.18);

        fluid_entity_collection_obj.step_2_set_mass_flowrate(
            mass_flowrate, 0);
        
        fluid_entity_collection_obj.step_2_set_mass_flowrate(
            mass_flowrate, 1);

        fluid_entity_collection_obj.step_2_set_mass_flowrate(
            mass_flowrate, 2);

        assert_eq!(0.180, 
                   fluid_entity_collection_obj.mass_flowrate_vec[0].value);
        assert_eq!(0.180, 
                   fluid_entity_collection_obj.mass_flowrate_vec[1].value);
        assert_eq!(0.180, 
                   fluid_entity_collection_obj.mass_flowrate_vec[2].value);

        use uom::si::power::watt;

        let pipe1_heat_rate = 
            Power::new::<watt>(100.0);
        let pipe2_heat_rate = 
            Power::new::<watt>(-20.0);
        let pipe3_heat_rate = 
            Power::new::<watt>(-80.0);

        let pipe_work_input_rate = 
            Power::new::<watt>(0.0);

        // let's set the heat and work

        fluid_entity_collection_obj.step_3_set_work_input(
            pipe_work_input_rate, 0);

        fluid_entity_collection_obj.step_3_set_work_input(
            pipe_work_input_rate, 1);

        fluid_entity_collection_obj.step_3_set_work_input(
            pipe_work_input_rate, 2);

        fluid_entity_collection_obj.step_4_set_heat_input(
            pipe1_heat_rate, 0);

        fluid_entity_collection_obj.step_4_set_heat_input(
            pipe2_heat_rate, 1);

        fluid_entity_collection_obj.step_4_set_heat_input(
            pipe3_heat_rate, 2);

        assert_eq!(100.0, 
                   fluid_entity_collection_obj.heat_input_vec[0].value);
        assert_eq!(-20.0, 
                   fluid_entity_collection_obj.heat_input_vec[1].value);
        assert_eq!(-80.0, 
                   fluid_entity_collection_obj.heat_input_vec[2].value);

        assert_eq!(0.0, 
                   fluid_entity_collection_obj.work_input_vec[0].value);
        assert_eq!(0.0, 
                   fluid_entity_collection_obj.work_input_vec[1].value);
        assert_eq!(0.0, 
                   fluid_entity_collection_obj.work_input_vec[2].value);
        // now let's do all the calculation and advance the timestep

        // current bug is here,
        // i get negative enthalpies,
        // already ruled out that setting heat and work input is wrong
        // however, the initial enthalpy of the system at setup
        // may be a cause with using 
        // export RUST_BACKTRACE=1 
        // as an environment variable

        fluid_entity_collection_obj.
            step_5_calculate_all_outlet_enthalpies_and_temperatures();
        fluid_entity_collection_obj.step_6_calculate_inlet_temperatures();
        fluid_entity_collection_obj.step_7_advance_timestep();


        // Now we need to assert that the outlet temperatures
        // are same as the new temperatures
        //
        // for reference, the outlet temperatures are as follows
        // T_new1 = 305.91 K
        // T_new2 = 298.8 K
        // T_new3 = 295.2 K


        let temp_1_val = 305.91_f64;
        let temp_2_val = 298.8_f64;
        let temp_3_val = 295.2_f64;

        pipe1 = fluid_entity_collection_obj.fluid_entity_vector[0].clone();
        pipe2 = fluid_entity_collection_obj.fluid_entity_vector[1].clone();
        pipe3 = fluid_entity_collection_obj.fluid_entity_vector[2].clone();

        approx::assert_relative_eq!(
            temp_1_val,
            pipe1.fluid_parameters.temperature_data.outlet_temp_old.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_2_val,
            pipe2.fluid_parameters.temperature_data.outlet_temp_old.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_3_val,
            pipe3.fluid_parameters.temperature_data.outlet_temp_old.value,
            max_relative=0.01);

        // likewise, let's assert the old inlet temperatures too

        approx::assert_relative_eq!(
            temp_1_val,
            pipe2.fluid_parameters.temperature_data.inlet_temp_old.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_2_val,
            pipe3.fluid_parameters.temperature_data.inlet_temp_old.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_3_val,
            pipe1.fluid_parameters.temperature_data.inlet_temp_old.value,
            max_relative=0.01);
        
    }

}


/// A structure to help calculate enthalpies and temperatures
/// for each time step. Experimental version 1
///
/// The fluid entity collection struct or class helps you to
/// calculate enthalpies and temperatures of Dowtherm A
/// pipes given a vector of work inputs, mass rate (flowrate)
/// inputs and heat inputs into the fluid control volume
///
///
/// Note: backwards flow is NOT supported yet.
///
/// The basic idea is to setup the collection first, 
/// this helps you set the initial temperatures,
/// add pipe objects, connect them in sequence
///
///
/// here's how you may use it:
///
/// There are three setup steps which you need to use:
/// ```rust
///
///
/// // start with some important imports
///
/// extern crate approx;
/// use heat_transfer_rust::ControlVolumeCalculations::Sandbox::
///     v2_IterativeHeatFluxTherminolPipe;
/// use heat_transfer_rust::ControlVolumeCalculations::ExplicitCalculations::
///     FluidEntityCollectionV1;
/// use heat_transfer_rust:: ControlVolumeCalculations::FluidEntity_StructsAndTraits::
///     FluidEntityInitialisationSteps;
/// use uom::si::f64::*;
/// use uom::si::time::second;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use uom::si::volume::cubic_meter;
///
///
///
///
/// // instantiate a new fluid entity collection object
///
/// let mut fluid_entity_collection_obj = 
///     FluidEntityCollectionV1::new();
///
///
///
/// // Setup Step 0: initiate timestep and global initial temp
///
/// let timestep = Time::new::<second>(0.1_f64);
/// let initial_global_temp = ThermodynamicTemperature::
///     new::<kelvin>(300_f64);
/// let fluid_volume = Volume::new::<cubic_meter>(
///     0.01_f64.powf(3_f64));
///
///
/// fluid_entity_collection_obj.setup_step_0_set_timestep_and_initial_temp(
///     timestep, 
///     initial_global_temp);
///
///
/// // Setup Step 1: we add 3 pipes
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe1".to_string(), 
///     fluid_volume);
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe2".to_string(), 
///     fluid_volume);
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe3".to_string(), 
///     fluid_volume);
///
/// 
/// // Setup Step 2: connect the pipes
/// // in this example, i connect like so:
///
///
/// // pipe1 --> pipe2 --> pipe3,
/// // pipe3 then connects back to pipe1
///
/// // to connect however, you need to know the index of the pipe
/// // the index is usually determined by the order you add
/// // the pipes to the fluid_entity_collection object
///
/// // pipe1 will be index 0,
/// // pipe2 will be index 1,
/// // pipe3 will be index 2
/// // and so on
///
///
/// // to check the index, it's kind of not user friendly,
/// // but you have to check the the list or vector of names
/// // get the index and check what the name of it is
/// // this user unfriendliness 
/// // can and should be sorted in future versions.
///
/// assert_eq!("pipe1".to_string(), 
/// fluid_entity_collection_obj.component_name_vec[0]);
/// 
/// assert_eq!("pipe2".to_string(), 
/// fluid_entity_collection_obj.component_name_vec[1]);
///
/// assert_eq!("pipe3".to_string(), 
/// fluid_entity_collection_obj.component_name_vec[2]);
///
/// // now that we know the index, we can connect the pipes by
/// // index
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(0, 1);
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(1, 2);
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(2, 0);
///
/// // now your setup is complete!
///
/// // end of setup demo
///
///
/// ```
///
///
/// Now that you have setup the problem, the next thing is to run it!
/// ```rust
///
///
/// // start with some important imports
///
/// extern crate approx;
/// use heat_transfer_rust::ControlVolumeCalculations::Sandbox::
///     v2_IterativeHeatFluxTherminolPipe;
/// use heat_transfer_rust::ControlVolumeCalculations::ExplicitCalculations::
///     FluidEntityCollectionV1;
/// use heat_transfer_rust:: ControlVolumeCalculations::FluidEntity_StructsAndTraits::
///     FluidEntityInitialisationSteps;
/// use uom::si::f64::*;
/// use uom::si::time::second;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use uom::si::volume::cubic_meter;
/// use uom::si::mass_rate::kilogram_per_second;
/// use uom::si::power::watt;
///
///
///
/// // instantiate a new fluid entity collection object
///
/// let mut fluid_entity_collection_obj = 
///     FluidEntityCollectionV1::new();
///
///
///
/// // Setup Step 0: initiate timestep and global initial temp
///
/// let timestep = Time::new::<second>(0.1_f64);
/// let initial_global_temp = ThermodynamicTemperature::
///     new::<kelvin>(300_f64);
/// let fluid_volume = Volume::new::<cubic_meter>(
///     0.01_f64.powf(3_f64));
///
///
/// fluid_entity_collection_obj.setup_step_0_set_timestep_and_initial_temp(
///     timestep, 
///     initial_global_temp);
///
///
/// // Setup Step 1: we add 3 pipes
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe1".to_string(), 
///     fluid_volume);
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe2".to_string(), 
///     fluid_volume);
///
/// fluid_entity_collection_obj.setup_step_1_add_new_component(
///     "pipe3".to_string(), 
///     fluid_volume);
///
/// 
/// // Setup Step 2: connect the pipes
/// // in this example, i connect like so:
///
///
/// // pipe1 --> pipe2 --> pipe3,
/// // pipe3 then connects back to pipe1
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(0, 1);
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(1, 2);
///
/// fluid_entity_collection_obj.
///     setup_step_2_connect_inlet_and_outlet_pipe(2, 0);
///
/// // now your setup is complete, we can move on to
/// // proper timestep calcs
///
/// // Let's now run through a sequence of calculation.
/// // Assume that for this timestep, 100W of heat is added
/// // to pipe1, 
/// // -80W of heat to pipe2
/// // and 
/// // -20W of heat to pipe3
/// // no work is done on the pipes
/// // and the mass flowrate through each pipe is
/// // 0.18 kg/s
///
/// // we can first initiate some of these
/// // power values first
///
/// let pipe1_heat_rate = 
///     Power::new::<watt>(100.0);
/// let pipe2_heat_rate = 
///     Power::new::<watt>(-20.0);
/// let pipe3_heat_rate = 
///     Power::new::<watt>(-80.0);
/// let pipe_work_input_rate = 
///     Power::new::<watt>(0.0);
/// 
///
/// let mass_flowrate = 
///    MassRate::new::<kilogram_per_second>(0.18);
/// 
/// fluid_entity_collection_obj.
///     step_1_calculate_current_timestep_temp_enthalpies();
///
/// fluid_entity_collection_obj.
/// step_2_set_mass_flowrate(mass_flowrate,0);
/// fluid_entity_collection_obj.
/// step_2_set_mass_flowrate(mass_flowrate,1);
/// fluid_entity_collection_obj.
/// step_2_set_mass_flowrate(mass_flowrate,2);
///
/// fluid_entity_collection_obj.
/// step_3_set_work_input(pipe_work_input_rate,0);
/// fluid_entity_collection_obj.
/// step_3_set_work_input(pipe_work_input_rate,1);
/// fluid_entity_collection_obj.
/// step_3_set_work_input(pipe_work_input_rate,2);
///
/// fluid_entity_collection_obj.
/// step_4_set_heat_input(pipe1_heat_rate,0);
/// fluid_entity_collection_obj.
/// step_4_set_heat_input(pipe2_heat_rate,1);
/// fluid_entity_collection_obj.
/// step_4_set_heat_input(pipe3_heat_rate,2);
///
/// fluid_entity_collection_obj.
/// step_5_calculate_all_outlet_enthalpies_and_temperatures();
/// fluid_entity_collection_obj.
/// step_6_calculate_inlet_temperatures();
/// fluid_entity_collection_obj.
/// step_7_advance_timestep();
///
///
/// // if everything works out right,
/// // the outlet temperatures of each pipe
///
/// let pipe1_outlet_temp_val = 305.91_f64;
/// let pipe2_outlet_temp_val = 298.8_f64;
/// let pipe3_outlet_temp_val = 295.2_f64;
///
/// let mut pipe1 = fluid_entity_collection_obj.fluid_entity_vector[0].clone();
/// let mut pipe2 = fluid_entity_collection_obj.fluid_entity_vector[1].clone();
/// let mut pipe3 = fluid_entity_collection_obj.fluid_entity_vector[2].clone();
///
/// approx::assert_relative_eq!(
///     pipe1_outlet_temp_val,
///     pipe1.fluid_parameters.temperature_data.outlet_temp_old.value,
///     max_relative=0.01);
///
/// approx::assert_relative_eq!(
///     pipe2_outlet_temp_val,
///     pipe2.fluid_parameters.temperature_data.outlet_temp_old.value,
///     max_relative=0.01);
///
/// approx::assert_relative_eq!(
///     pipe3_outlet_temp_val,
///     pipe3.fluid_parameters.temperature_data.outlet_temp_old.value,
///     max_relative=0.01);
///
/// // additionally, pipe3's outlet is the same as pipe1's inlet
/// // pipe1's outlet is the same as pipe2's inlet
/// // pipe2's outlet is the same as pipe3's inlet
/// // this was based on the way we connected our pipes
///
/// approx::assert_relative_eq!(
///     pipe1_outlet_temp_val,
///     pipe2.fluid_parameters.temperature_data.inlet_temp_old.value,
///     max_relative=0.01);
///
/// approx::assert_relative_eq!(
///     pipe2_outlet_temp_val,
///     pipe3.fluid_parameters.temperature_data.inlet_temp_old.value,
///     max_relative=0.01);
///
/// approx::assert_relative_eq!(
///     pipe3_outlet_temp_val,
///     pipe1.fluid_parameters.temperature_data.inlet_temp_old.value,
///     max_relative=0.01);
///
/// ```
///
///
///
/// to update in future version:
///
/// (1) creating pipes with strings should be easier, so that we do not have
/// to use the .to_string() method
///
/// (2) connecting pipes by name should be easier, so we don't have to use
/// pipe indices, but rather pipe names. So
/// we should be able to use both names and indices to connect the pipes
///
/// (3) backwards flow, bidirectional flow
///
/// (4) iter_mut() methods in for loops for parallel
///
/// (5) the cooling rate of the pipe should be automatically determined,
/// doing it manually every timestep is quite tedious
///
/// (6) point 5 same goes for work input, and mass flowrate, the user
/// should specify if there is a standard heat input rate at the heater
/// and the rest should be autocalculated
///
/// (7) generic fluid entity class that is able to perform all these
/// methods via traits no matter the underlying type (good for heat
/// exchangers and other components that have different nusselt number
/// configurations, or connected to some underlying heat structure)
///
/// (8) for traits with default implementations, we cannot act upon
/// the struct directly, however we can interact with the struct
/// via get/set methods which are to be implemented. Will probably
/// want a trait with get/set methods for interfacing with variables
/// and then generic methods which can then use those get and set
/// methods.
///
///
///
///
///
pub struct FluidEntityCollectionV1 {

    pub current_max_index: usize,
    pub fluid_entity_vector: Vec<v2_IterativeHeatFluxTherminolPipe>,
    pub inlet_temp_vec: Vec<ThermodynamicTemperature>,
    pub outlet_temp_vec: Vec<ThermodynamicTemperature>,

    pub heat_input_vec: Vec<Power>,
    pub work_input_vec: Vec<Power>,
    pub mass_flowrate_vec: Vec<MassRate>,

    pub component_name_vec: Vec<String>,

    pub timestep: Time,
    pub initial_global_temp: ThermodynamicTemperature,

}

use crate::ControlVolumeCalculations::FluidEntity_StructsAndTraits::*;

impl FluidEntityCollectionV1 {

    /// default constructor
    /// sets timestep at 0.1s by default
    pub fn new() -> Self {

        return Self { 
            current_max_index: 0, 
            fluid_entity_vector: vec![], 
            inlet_temp_vec: vec![], 
            outlet_temp_vec: vec![], 
            heat_input_vec: vec![], 
            work_input_vec: vec![], 
            mass_flowrate_vec: vec![], 
            component_name_vec: vec![],

            /// Default timestep is 0.1s
            timestep: Time::new::<second>(0.1),
            /// Default initial global temp is 300k or about 27C
            initial_global_temp: ThermodynamicTemperature::new::
                <kelvin>(300.0),
        }

    }

    pub fn setup_step_0_set_timestep_and_initial_temp(
        &mut self,
        timestep: Time,
        initial_global_temp: ThermodynamicTemperature) {

        self.timestep = timestep;
        self.initial_global_temp = initial_global_temp;

        
    }



    /// Adds a new fluid component or fluid entity
    pub fn setup_step_1_add_new_component(
        &mut self,
        name: String,
        fluid_volume: Volume
        ) {

        // first let me push the name of the fluid entity up
        self.component_name_vec.push(name);

        let mut new_pipe = v2_IterativeHeatFluxTherminolPipe::new();

        let fluid_entity_index = self.current_max_index;

        // we make a new fluid entity
        new_pipe.step_0_set_timestep_and_initial_temperatures(
            self.timestep,
            self.initial_global_temp,
            fluid_volume,
            fluid_entity_index);

        // i will push the default new pipe with fluid volume
        // and push the default initial temperatures into the
        // temperature inlet and outlet vectors
        self.fluid_entity_vector.push(new_pipe);
        self.inlet_temp_vec.push(self.initial_global_temp.clone());
        self.outlet_temp_vec.push(self.initial_global_temp.clone());

        // i will then push the default work done and heat input into
        // these vectors

        self.work_input_vec.push(
            Power::new::<watt>(0.0));
        self.heat_input_vec.push(
            Power::new::<watt>(0.0));


        // the mass flowrate is also by default set to zero

        self.mass_flowrate_vec.push(
            MassRate::new::<kilogram_per_second>(0.0));


        // add 1 to the maximum current index

        self.current_max_index = self.current_max_index + 1;

        return;

    }

    /// setup step : connect inlet and outlet of pipe
    pub fn setup_step_2_connect_inlet_and_outlet_pipe(
        &mut self,
        connect_to_pipe_outlet_index: usize,
        connect_to_pipe_inlet_index: usize){


        // Basically in this function, i cannot use borrow
        // two mutable versions of the component vector and then
        // change values in them
        // i have to make a copy of the front and back pipe
        // and then use those to perform the value changes
        // the other way of course, is to change the value indices manually


        let mut pipe_back = 
            self.fluid_entity_vector[connect_to_pipe_outlet_index].clone();

        let mut pipe_front = 
            self.fluid_entity_vector[connect_to_pipe_inlet_index].clone();

        pipe_front.step_1_connect_to_component_inlet(
            &mut self.fluid_entity_vector[connect_to_pipe_outlet_index]);

        pipe_back.step_2_conenct_to_component_outlet(
            &mut self.fluid_entity_vector[connect_to_pipe_inlet_index]);

    }

    /// Now we are going into running the simulation
    /// the first step is to calculate the current timestep 
    /// temperature and enthalpies

    pub fn step_1_calculate_current_timestep_temp_enthalpies(
        &mut self) {

        // start the for loop
        let max_vec_index_plus_one = 
            self.fluid_entity_vector.len();

        for i in 0..max_vec_index_plus_one {
            self.fluid_entity_vector[i].
                step_1_calculate_current_timestep_temp_enthalpies();
        }
    }

    /// Step 2: set mass flowrate for a component with 
    /// index i
    pub fn step_2_set_mass_flowrate(
        &mut self,
        mass_flowrate: MassRate,
        component_index: usize){

        self.mass_flowrate_vec[component_index] = mass_flowrate.clone();
    }

    /// Step 3: set work input vector for component with 
    /// index i
    pub fn step_3_set_work_input(
        &mut self,
        work_input: Power,
        component_index: usize){

        self.work_input_vec[component_index] = work_input.clone();
    }

    /// Step 4: set heat input vector for component with 
    /// index i
    pub fn step_4_set_heat_input(
        &mut self,
        heat_input: Power,
        component_index: usize){

        self.heat_input_vec[component_index] = heat_input.clone();
    }

    /// Step 5: calculate outlet enthalpy
    /// in a serial manner (not worrying about parallel
    /// computation with rayon yet)
    pub fn step_5_calculate_all_outlet_enthalpies_and_temperatures(
        &mut self) {

        // start the for loop
        let max_vec_index_plus_one = 
            self.fluid_entity_vector.len();

        for i in 0..max_vec_index_plus_one {

            let heat_input_into_fluid : Power = 
                self.heat_input_vec[i];

            let work_done_on_fluid : Power = 
                self.work_input_vec[i];

            let mass_flowrate: MassRate = 
                self.mass_flowrate_vec[i];

            self.fluid_entity_vector[i].
                step_2_calculate_new_outlet_enthalpy(
                    heat_input_into_fluid,
                    work_done_on_fluid,
                    self.timestep,
                    mass_flowrate);

            self.fluid_entity_vector[i].
                step_3_calculate_new_outlet_temperature();
        }
        return;

    }

    /// Step 6: calculate inlet temperatures and assign them
    /// to appropriate vectors
    pub fn step_6_calculate_inlet_temperatures(
        &mut self){

        // first let's clear up the inlet and outlet temp vector

        self.inlet_temp_vec.clear();
        self.outlet_temp_vec.clear();
        // now let's update all the outlet temperatures

        let max_vec_index_plus_one = 
            self.fluid_entity_vector.len();

        // first we update the outlet temp vector

        for i in 0..max_vec_index_plus_one {

            // first let's obtain the outlet temperature
            //
            let fluid_component =
                self.fluid_entity_vector[i].clone();

            let fluid_component_outlet_temperature =
                fluid_component.fluid_parameters.
                temperature_data.outlet_temp_new;

            // we'll introduce it into the vector

            self.outlet_temp_vec.push(fluid_component_outlet_temperature);

            // of course, we can set the outlet temperature here outright,
            // but i'll leave it for later


        }


        // now we update the inlet temp vector
        //
        for i in 0.. max_vec_index_plus_one {

            let fluid_component =
                self.fluid_entity_vector[i].clone();

            let fluid_component_inlet_index: usize = 
                fluid_component.fluid_parameters.
                index_data.inlet_fluid_entity_index;

            // second, we get the outlet temperature of the fluid
            // component connected to the back of this fluid
            // component
            //
            // This is actually the inlet temperature

            let fluid_component_inlet_temperature =
                self.outlet_temp_vec[fluid_component_inlet_index];

            // third, let's push this temperature to
            // the inlet temperature vector

            self.inlet_temp_vec.push(fluid_component_inlet_temperature);

        }


        // now we update the inlet and outlet temperatures
        // of each object
        //

        for i in 0..max_vec_index_plus_one {

            self.fluid_entity_vector[i].fluid_parameters.
                temperature_data.inlet_temp_new 
                = self.inlet_temp_vec[i].clone();

            self.fluid_entity_vector[i].fluid_parameters.
                temperature_data.outlet_temp_new =
                self.outlet_temp_vec[i].clone();
        }

        // and we are done!
        return;

    }


    /// Step 7: Advance timestep
    /// means i set the old temperature values to that of the new
    /// temperature
    pub fn step_7_advance_timestep(
        &mut self){

        let max_vec_index_plus_one = 
            self.fluid_entity_vector.len();

        for i in 0..max_vec_index_plus_one {
            self.fluid_entity_vector[i].
                step_6_update_current_timestep_temperatures();
        }
    }


}

/// For explicit calculations in general
///
/// we use the equation:
///
/// h_new = h_old + deltaT * (H_in - H_out +
/// Q + W)
///
/// What we need here:
///
/// we take in the current T_in and T_out of
/// the pipe or thermal component
///
/// old and new as well
///
/// and T_sys old and new
///
/// (1) To get enthalpy, I will then need to convert
/// these temperatures to enthalpy, so i will need methods
/// for that
///
///
/// This trait helps the developer run through the steps
/// of enthalpy calculation
///
///
pub trait v2_ExplicitCalculationSteps {


    /// First Step: calculate enthalpies and bulk fluid temp
    /// from temperatures
    fn step_1_calculate_current_timestep_temp_enthalpies(
        &mut self);

    /// Second Step: calculate new outlet enthalpy from available
    /// enthalpies, heat loss/gain and work done rates
    /// will probably require timestep also
    fn step_2_calculate_new_outlet_enthalpy(
        &mut self, 
        heat_supplied_to_fluid: Power,
        work_done_on_fluid: Power,
        timestep: Time,
        fluid_mass_flowrate: MassRate);

    /// third step: calculate new outlet
    /// temperature based on new outlet enthalpy
    ///
    /// we will then obtain the inlet and outlet
    /// temperatures based on the outlet temperatures
    /// of each pipe using some matrix solver
    /// which should give us a vector of inlet
    /// temperatures
    fn step_3_calculate_new_outlet_temperature(
        &mut self) -> ThermodynamicTemperature;

    /// after finding the vector of inlet temperatures
    /// we should either give the thermodynamic temperautre
    /// and put it into the function
    /// 
    /// or make a copy of the vector of inlet temperature
    /// and feed it to this struct so that it can edit 
    /// its own data
    fn step_4_set_inlet_temperature(
        &mut self,
        new_inlet_temperature: ThermodynamicTemperature);

    /// after finding the vector of inlet temperatures
    /// we can next map the inlet temperatures to the
    /// proper outlet temperature vector and also feed it
    /// in
    fn step_5_set_outlet_temperature(
        &mut self,
        new_outlet_temperature: ThermodynamicTemperature);


    /// now that we have all the required information,
    /// we can set the old temperatures to the values of
    /// the new temperatures
    fn step_6_update_current_timestep_temperatures(
        &mut self);
}


/// old version of explicit calculation steps,
/// i thought it was okay to calculate new system temperature
/// T_sys = (T_in+T_out)/2
///
/// it turns out that this method does not conserver energy
/// and so i will have to use the well mixed assumption
pub trait v1_ExplicitCalculationSteps {


    /// First Step: calculate enthalpies and bulk fluid temp
    /// from temperatures
    fn step_1_calculate_current_timestep_temp_enthalpies(
        &mut self);

    /// Second Step: calculate new system enthalpy from available
    /// enthalpies, heat loss/gain and work done rates
    /// will probably require timestep also
    fn step_2_calculate_new_system_enthalpy(
        &mut self, 
        heat_supplied_to_fluid: Power,
        work_done_on_fluid: Power,
        timestep: Time,
        fluid_mass_flowrate: MassRate);

    /// third step: calculate new system
    /// temperature based on new system enthalpy
    ///
    /// we will then obtain the inlet and outlet
    /// temperatures based on the system temperatures
    /// of each pipe using some matrix solver
    /// which should give us a vector of inlet
    /// temperatures
    fn step_3_calculate_new_system_temperature(
        &mut self) -> ThermodynamicTemperature;

    /// after finding the vector of inlet temperatures
    /// we should either give the thermodynamic temperautre
    /// and put it into the function
    /// 
    /// or make a copy of the vector of inlet temperature
    /// and feed it to this struct so that it can edit 
    /// its own data
    fn step_4_set_inlet_temperature(
        &mut self,
        new_inlet_temperature: ThermodynamicTemperature);

    /// after finding the vector of inlet temperatures
    /// we can next map the inlet temperatures to the
    /// proper outlet temperature vector and also feed it
    /// in
    fn step_5_set_outlet_temperature(
        &mut self,
        new_outlet_temperature: ThermodynamicTemperature);


    /// now that we have all the required information,
    /// we can set the old temperatures to the values of
    /// the new temperatures
    fn step_6_update_current_timestep_temperatures(
        &mut self);
}


