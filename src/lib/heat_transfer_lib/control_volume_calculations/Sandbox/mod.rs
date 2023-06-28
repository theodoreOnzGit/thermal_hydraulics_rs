#[warn(missing_docs)]

#[cfg(test)]
mod sandbox_therminol_dowtherm_pipes {
    /// from this sandbox, i am experimenting with 
    /// using a struct and its methods to build my therminolPipes
    /// and run through their calculation steps at every timestep
    ///
    /// so far i've learnt that:
    ///
    /// (1) it's better to calculate the outlet temperature directly
    /// rather than go for fluid temperature; the using
    /// T_fluid = (T_in+T_out)/2 does not conserve energy
    /// --> (done and tested 19 dec 2022)
    ///
    /// 
    ///
    /// (2) let the mutable inlet temperature and outlet temperautre
    /// vectors be a a part of the struct, same goes for the calculator
    /// --> (done and tested 19 dec 2022)
    ///
    ///
    /// (3) index numbers and vectors to keep track of inlet and outlet
    /// temperatures seems pretty okay! --> (no change needed)
    ///
    ///
    ///
    /// (4) setting inlet and outlet temperatures should be the same
    /// step 
    /// --> (done and tested 19 dec 2022)
    ///
    ///
    ///
    /// (5) use usize instead of i32 to index vectors, though i32 
    /// can be readily converted into usize. --> (changed 1128 hrs 16 Dec 2022)
    ///
    ///
    ///
    /// (6) power, massfow and work done vectors can be initiated in the struct
    /// as well
    /// --> (done and tested 19 dec 2022)
    ///
    ///
    /// 
    /// (7) we can set a global timestep and global time in our 
    /// struct
    /// --> (done and tested 19 dec 2022)
    ///
    ///
    /// (8) we may want to have get and set methods for power/work done
    /// and mass flow
    /// --> (set methods only 19 dec 2022 done and tested)
    ///
    ///
    ///
    /// for future consideration:
    ///
    /// (1) the inputs for the calculator function may need to extend
    /// beyond simple pipes if we are to include heat exchangers
    /// and other things
    ///
    /// (2) we might want to implement an ambient .cool() function
    /// in order to represent heat loss to environment
    ///
    /// (3) bidirectional flow...
    #[test]
    pub fn sandbox_autobuild_conenct_and_autocalculation() {
        extern crate approx;

        use crate::ControlVolumeCalculations::Sandbox::
            v1_IterativeHeatFluxTherminolPipe;

        use crate:: ControlVolumeCalculations::FluidEntity_StructsAndTraits::
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

        // first thing first, i want to streamline the pipe creation process
        // (and of course, eventually place those hydrodynamic parameters as well,
        // but will probably think about that later)
        // 
        // and i also want to use a for loop to do these calculation steps
        // hence, putting these in a vector form will be good
        // so i will have 3 vectors,
        // one vector for inlet temperatures
        // one vector for outlet temperatures
        // one vector for all the v1_IterativeHeatFluxTherminolPipe

        let mut therminolPipeVec: Vec<v1_IterativeHeatFluxTherminolPipe> 
            = vec![];

        let mut inlet_temp_vec: Vec<ThermodynamicTemperature>
            = vec![];

        let mut outlet_temp_vec: Vec<ThermodynamicTemperature>
            = vec![];

        pub struct HeatFluxPipeFactory {
            pub current_max_index: usize,
        }

        impl HeatFluxPipeFactory {

            pub fn new() -> Self {

                return Self { current_max_index: 0 };

            }

            pub fn add_new_component(
                &mut self,
                timestep: Time,
                initial_global_temp: ThermodynamicTemperature,
                fluid_volume: Volume,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>,
                inlet_temp_vec: &mut Vec<ThermodynamicTemperature>,
                outlet_temp_vec: &mut Vec<ThermodynamicTemperature>
                ) {

                let mut new_pipe = v1_IterativeHeatFluxTherminolPipe::new();

                let fluid_entity_index = self.current_max_index;

                new_pipe.step_0_set_timestep_and_initial_temperatures(
                    timestep,
                    initial_global_temp,
                    fluid_volume,
                    fluid_entity_index);

                generic_component_vec.push(new_pipe);
                inlet_temp_vec.push(initial_global_temp);
                outlet_temp_vec.push(initial_global_temp);

                self.current_max_index = fluid_entity_index + 1;

            }

        }

        let mut factory_obj = HeatFluxPipeFactory::new();

        // so this is pipe 1

        factory_obj.add_new_component(
            timestep,
            initial_global_temp,
            fluid_volume,
            &mut therminolPipeVec,
            &mut inlet_temp_vec,
            &mut outlet_temp_vec);

        // only issue here, no name for each component, may be difficult
        // to differentiate


        // so this is pipe 2

        factory_obj.add_new_component(
            timestep,
            initial_global_temp,
            fluid_volume,
            &mut therminolPipeVec,
            &mut inlet_temp_vec,
            &mut outlet_temp_vec);


        // so this is pipe 3

        factory_obj.add_new_component(
            timestep,
            initial_global_temp,
            fluid_volume,
            &mut therminolPipeVec,
            &mut inlet_temp_vec,
            &mut outlet_temp_vec);

        // the next step is to connect each pipe
        // so i can try giving a factory method which helps to connect
        // two pipe objects
        // for this i will need a pipe vector, 

        // the index of the pipe which
        // we are interested to connect at the outlet
        // and the pipe we are interested to connect at the inlet

        impl HeatFluxPipeFactory {

            pub fn connect_inlet_and_outlet_pipe(
                &self,
                connect_to_pipe_outlet_index: usize,
                connect_to_pipe_inlet_index: usize,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>){


                // Basically in this function, i cannot use borrow
                // two mutable versions of the component vector and then
                // change values in them
                // i have to make a copy of the front and back pipe
                // and then use those to perform the value changes
                // the other way of course, is to change the value indices manually


                let mut pipe_back = 
                    generic_component_vec[connect_to_pipe_outlet_index].clone();

                let mut pipe_front = 
                    generic_component_vec[connect_to_pipe_inlet_index].clone();

                pipe_front.step_1_connect_to_component_inlet(
                    &mut generic_component_vec[connect_to_pipe_outlet_index]);

                pipe_back.step_2_conenct_to_component_outlet(
                    &mut generic_component_vec[connect_to_pipe_inlet_index]);

            }


        }

        // I can use the above method to connect my pipes!
        // i have it this way:
        // 1 -> 2 -> 3 
        // and 3 connects back to 1 in a circular fashion

        // however pipe 1 has an index of 0,
        // pipe 2 has an index of 1 and so on...

        factory_obj.connect_inlet_and_outlet_pipe(0,1,
                                                  &mut therminolPipeVec);

        factory_obj.connect_inlet_and_outlet_pipe(1,2,
                                                  &mut therminolPipeVec);

        factory_obj.connect_inlet_and_outlet_pipe(2,0,
                                                  &mut therminolPipeVec);


        // in doing so, i have connected all the pipes and these
        // pipes will know the correct index

        // all right, now let's take pipe1 for example,
        // its own index should be 0
        // the pipe connected to its back is index 2,
        // the pipe connected to its front is index 1
        // i'm going to repeat it all for the pipes involved

        let pipe1 = therminolPipeVec[0].clone();

        assert_eq!(0, pipe1.fluid_parameters.index_data.fluid_entity_index);

        assert_eq!(2, pipe1.fluid_parameters.index_data.
                   inlet_fluid_entity_index);

        assert_eq!(1, pipe1.fluid_parameters.index_data.
                   outlet_fluid_entity_index);

        let pipe2 = therminolPipeVec[1].clone();

        assert_eq!(1, pipe2.fluid_parameters.index_data.
                   fluid_entity_index);

        assert_eq!(0, pipe2.fluid_parameters.index_data.
                   inlet_fluid_entity_index);

        assert_eq!(2, pipe2.fluid_parameters.index_data.
                   outlet_fluid_entity_index);

        let pipe3 = therminolPipeVec[2].clone();

        assert_eq!(2, pipe3.fluid_parameters.index_data.
                   fluid_entity_index);
        assert_eq!(1, pipe3.fluid_parameters.index_data.
                   inlet_fluid_entity_index);
        assert_eq!(0, pipe3.fluid_parameters.index_data.
                   outlet_fluid_entity_index);


        // now that we have connected the pipes, we can start
        // calculating
        // we can perhaps create two vectors,
        // one for work done, one for heat input in fluid
        // one for work done on fluid
        //
        use crate::ControlVolumeCalculations::ExplicitCalculations::
            v1_ExplicitCalculationSteps;

        impl HeatFluxPipeFactory {

            pub fn step_1_calculate_current_timestep_temp_enthalpies(&self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>){

                // start function

                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {
                    generic_component_vec[i].
                        step_1_calculate_current_timestep_temp_enthalpies();

                }


                // end of function
            }

            pub fn step_2_calculate_new_system_enthalpy(
                &self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>,
                heat_input_vec: Vec<Power>,
                work_input_vec: Vec<Power>,
                mass_flowrate_vec: Vec<MassRate>,
                timestep: Time) {

                // start function
                
                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {

                    // first let's get the heat rate:
                    //
                    let heat_input_into_fluid : Power = 
                        heat_input_vec[i];

                    // now work done rate:

                    let work_input_into_fluid : Power = 
                        work_input_vec[i];

                    // now mass flowrate:
                    //

                    let mass_flowrate: MassRate = 
                        mass_flowrate_vec[i];




                    generic_component_vec[i].
                        step_2_calculate_new_system_enthalpy(
                            heat_input_into_fluid,
                            work_input_into_fluid,
                            timestep,
                            mass_flowrate);

                }

                // end of function
            }

            pub fn step_3_calculate_new_system_temperature(
                &mut self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>){


                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {
                    generic_component_vec[i].
                        step_3_calculate_new_system_temperature();

                }



            }

        }

        // we can then use this function to help us with step 1
        //

        factory_obj.step_1_calculate_current_timestep_temp_enthalpies(
            &mut therminolPipeVec);

        // now for calculation
        let timestep : Time = pipe1.fluid_parameters.timestep;

        use uom::si::mass_rate::kilogram_per_second;

        let mass_flowrate: MassRate = 
            MassRate::new::<kilogram_per_second>(0.18);

        use uom::si::power::watt;
        
        let pipe_1_heat = Power::new::<watt>(100_f64);
        let pipe_2_heat = Power::new::<watt>(-20_f64);
        let pipe_3_heat = Power::new::<watt>(-80_f64);
        
        let work_done_on_pipe_rate = Power::new::<watt>(0_f64);

        // let's construct the vectors
        //

        let mut heat_input_vec : Vec<Power> =
            vec![];

        let mut work_input_vec : Vec<Power> =
            vec![];

        let mut mass_flowrate_vec: Vec<MassRate> = 
            vec![];

        // push values into vector

        heat_input_vec.push(pipe_1_heat);
        heat_input_vec.push(pipe_2_heat);
        heat_input_vec.push(pipe_3_heat);

        work_input_vec.push(work_done_on_pipe_rate);
        work_input_vec.push(work_done_on_pipe_rate);
        work_input_vec.push(work_done_on_pipe_rate);

        mass_flowrate_vec.push(mass_flowrate);
        mass_flowrate_vec.push(mass_flowrate);
        mass_flowrate_vec.push(mass_flowrate);


        // next thing, complete step 2

        factory_obj.step_2_calculate_new_system_enthalpy(
            &mut therminolPipeVec,
            heat_input_vec,
            work_input_vec,
            mass_flowrate_vec,
            timestep);

        // step 3

        factory_obj.step_3_calculate_new_system_temperature(
            &mut therminolPipeVec);


        // now at this point, there is still some correction to 
        // make with the v1_IterativeHeatFluxTherminolPipe
        //
        // since the enthalpy balances do not corroborate well
        // 
        // but the new fluid temperatures should now reflect
        // the outlet temperatures of the pipes
        // Hence, i will take the system temperatures of the 
        // pipe and assign it appropriately to the vector of
        // the outlet temperatures

        impl HeatFluxPipeFactory {

            pub fn step_4_set_inlet_temperature(
                &mut self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>,
                inlet_temp_vec: &mut Vec<ThermodynamicTemperature>,
                outlet_temp_vec: &mut Vec<ThermodynamicTemperature>){

                // first let's clear up the inlet and outlet temp vector

                inlet_temp_vec.clear();
                outlet_temp_vec.clear();

                // now let's update all the outlet temperatures

                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {
                    
                    // first let's obtain the outlet temperature
                    //
                    let fluid_component =
                        generic_component_vec[i].clone();

                    let fluid_component_outlet_temperature =
                        fluid_component.fluid_parameters.
                        temperature_data.fluid_temp_new;

                    // we'll introduce it into the vector

                    outlet_temp_vec.push(fluid_component_outlet_temperature);

                    // of course, we can set the outlet temperature here outright,
                    // but i'll leave it for later

                }

                // now that we've set the outlet temperature
                // vectors, we can use the indexing within each object
                // to assign them to the inlet temperatures

                for i in 0..=max_vec_index {

                    // first let's obtain the index
                    // of the component connected to this component's
                    // inlet
                    
                    let fluid_component =
                        generic_component_vec[i].clone();

                    let fluid_component_inlet_index : usize =
                        fluid_component.fluid_parameters.
                        index_data.inlet_fluid_entity_index
                        .try_into().unwrap();

                    // second, we get the outlet temperature of the fluid
                    // component connected to the back of this fluid
                    // component
                    //
                    // This is actually the inlet temperature

                    let fluid_component_inlet_temperature =
                        outlet_temp_vec[fluid_component_inlet_index];

                    // third, let's push this temperature to
                    // the inlet temperature vector

                    inlet_temp_vec.push(fluid_component_inlet_temperature);

                }

                // now we have set both inlet and outlet temperature vectors,
                // we can start assigning inlet temperatures

                for i in 0..=max_vec_index {

                    generic_component_vec[i].fluid_parameters.
                        temperature_data.inlet_temp_new 
                        = inlet_temp_vec[i].clone();

                }

                return;


            }

            pub fn step_5_set_outlet_temperature(
                &mut self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>,
                outlet_temp_vec: &mut Vec<ThermodynamicTemperature>){

                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {

                    generic_component_vec[i].fluid_parameters.
                        temperature_data.outlet_temp_new =
                        outlet_temp_vec[i].clone();

                }

                return;
            }

            // end of impl
        }


        factory_obj.step_4_set_inlet_temperature(
            &mut therminolPipeVec,
            &mut inlet_temp_vec,
            &mut outlet_temp_vec);

        factory_obj.step_5_set_outlet_temperature(
            &mut therminolPipeVec,
            &mut outlet_temp_vec);

        // now let's do some asserts

        let mut pipe1 = therminolPipeVec[0].clone();
        let mut pipe2 = therminolPipeVec[1].clone();
        let mut pipe3 = therminolPipeVec[2].clone();

        // for reference, the outlet temperatures are as follows
        // T_new1 = 305.91 K
        // T_new2 = 298.8 K
        // T_new3 = 295.2 K


        let temp_1_val = 305.91_f64;
        let temp_2_val = 298.8_f64;
        let temp_3_val = 295.2_f64;

        approx::assert_relative_eq!(
            temp_1_val,
            pipe1.fluid_parameters.temperature_data.outlet_temp_new.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_2_val,
            pipe2.fluid_parameters.temperature_data.outlet_temp_new.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_3_val,
            pipe3.fluid_parameters.temperature_data.outlet_temp_new.value,
            max_relative=0.01);

        // likewise, let's assert the new inlet temperatures too

        approx::assert_relative_eq!(
            temp_1_val,
            pipe2.fluid_parameters.temperature_data.inlet_temp_new.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_2_val,
            pipe3.fluid_parameters.temperature_data.inlet_temp_new.value,
            max_relative=0.01);

        approx::assert_relative_eq!(
            temp_3_val,
            pipe1.fluid_parameters.temperature_data.inlet_temp_new.value,
            max_relative=0.01);

        impl HeatFluxPipeFactory {
            
            pub fn step_6_update_current_timestep_temperatures(
                &mut self,
                generic_component_vec: &mut Vec<v1_IterativeHeatFluxTherminolPipe>){


                let max_vec_index = 
                    generic_component_vec.len() - 1;

                for i in 0..=max_vec_index {
                    generic_component_vec[i].
                        step_6_update_current_timestep_temperatures();

                }

                return;
            }

            // end of impl
        }

        // now we change the old temperature data (from last timestep)
        // to the new temperature data (current timestep)
        factory_obj.step_6_update_current_timestep_temperatures(
            &mut therminolPipeVec);

        pipe1 = therminolPipeVec[0].clone();
        pipe2 = therminolPipeVec[1].clone();
        pipe3 = therminolPipeVec[2].clone();

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

extern crate uom;
use uom::si::available_energy::joule_per_kilogram;
use uom::si::f64::*;
use uom::si::*;
use uom::si::heat_capacity::joule_per_kelvin;
use crate::ControlVolumeCalculations::ExplicitCalculations::*;
use crate::ControlVolumeCalculations::
FluidEntity_StructsAndTraits::*;


// here are imports for units
use uom::si::time::second;
use uom::si::volume::cubic_meter;
use uom::si::thermodynamic_temperature::kelvin;

extern crate fluid_mechanics_rust;
use fluid_mechanics_rust::therminol_component::
FluidProperties;
use fluid_mechanics_rust::therminol_component::
dowtherm_a_properties;


#[derive(Clone)]
pub struct v2_IterativeHeatFluxTherminolPipe {
    pub fluid_parameters: FluidEntityThermophysicalData,
}

impl v2_IterativeHeatFluxTherminolPipe {

    /// Constructor which creates the structure 
    /// ```rust
    /// 
    /// extern crate approx;
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// TherminolDowthermPipes::*;
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// Sandbox::*;
    ///
    /// use uom::si::f64::*;
    /// use uom::si::time::second;
    /// use uom::si::thermodynamic_temperature::kelvin;
    /// use uom::si::volume::cubic_meter;
    ///
    /// // first let's initialise a pipe using a timestep
    /// // global thermodynamic temperature and fluid volume
    /// 
    /// let timestep = Time::new::<second>(0.1_f64);
    /// let initial_global_temp = ThermodynamicTemperature::
    /// new::<kelvin>(300_f64);
    ///
    /// let fluid_volume = Volume::new::<cubic_meter>(
    /// 0.01_f64.powf(3_f64));
    ///
    /// let fluid_entity_index: usize = 4;
    ///
    /// // we are now going to initialise stuff
    ///
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// FluidEntity_StructsAndTraits::FluidEntityInitialisationSteps;
    ///
    /// let mut pipe1 = v2_IterativeHeatFluxTherminolPipe::new();
    ///
    /// pipe1.step_0_set_timestep_and_initial_temperatures(
    /// timestep,
    /// initial_global_temp,
    /// fluid_volume,
    /// fluid_entity_index);
    ///
    ///
    ///
    /// ```
    pub fn new() -> Self {

        // we'll need to make a few data structures first
        // with some default values


        let default_timestep: Time = 
            Time::new::<second>(0.1);

        let default_fluid_volume: Volume = 
            Volume::new::<cubic_meter>(0.02_f64.powf(3.0));

        // let's populate default index data
        
        let default_index: usize = 0;

        let default_index_data = 
            FluidEntityIndexData { 
                fluid_entity_index : default_index ,
                inlet_fluid_entity_index : default_index,
                outlet_fluid_entity_index : default_index,
            };


        // let's populate default temperature data
        //

        let default_temperature : ThermodynamicTemperature 
            = ThermodynamicTemperature::new::
            <kelvin>(310_f64);

        let default_temperature_data : PipeFluidTemperatureData  
            = PipeFluidTemperatureData { 
                inlet_temp_old: default_temperature, 
                inlet_temp_new: default_temperature, 
                outlet_temp_old: default_temperature, 
                outlet_temp_new: default_temperature, 
                fluid_temp_old: default_temperature, 
                fluid_temp_new: default_temperature,
            };

        // let's finally populate the enthalpy data

        let default_enthalpy : AvailableEnergy =
            AvailableEnergy::new::<joule_per_kilogram>(0_f64);

        let default_enthalpy_data : PipeFluidEnthalpyData =
            PipeFluidEnthalpyData { 
                inlet_enthalpy_old: default_enthalpy, 
                inlet_enthalpy_new: default_enthalpy, 
                outlet_enthalpy_old: default_enthalpy, 
                outlet_enthalpy_new: default_enthalpy, 
                fluid_enthalpy_old: default_enthalpy, 
                fluid_enthalpy_new: default_enthalpy 
            };

        let default_thermophysical_data : FluidEntityThermophysicalData
            = FluidEntityThermophysicalData { 
                index_data: default_index_data, 
                temperature_data: default_temperature_data, 
                enthalpy_data: default_enthalpy_data, 
                timestep: default_timestep, 
                fluid_volume: default_fluid_volume ,
            };

        return Self { 
            fluid_parameters: default_thermophysical_data 
        };

    }
}


impl FluidEntityInitialisationSteps for v2_IterativeHeatFluxTherminolPipe {
    /// step zero, here we set our initial conditions

    fn step_0_set_timestep_and_initial_temperatures(
        &mut self,
        timestep: Time,
        initial_global_temp: ThermodynamicTemperature,
        fluid_volume: Volume,
        fluid_entity_index: usize) -> Self {

        FluidEntityThermophysicalData::
            step_0_set_timestep_and_initial_temperatures(
                &mut self.fluid_parameters, 
                timestep, 
                initial_global_temp, 
                fluid_volume, 
                fluid_entity_index
                );

        return self.clone();
    }
    

    fn step_1_connect_to_component_inlet(
        &mut self,
        other_fluid_entity: &mut Self) {

        FluidEntityThermophysicalData::
            step_1_connect_to_component_inlet(
                &mut self.fluid_parameters, 
                &mut other_fluid_entity.fluid_parameters);
    }

    fn step_2_conenct_to_component_outlet(
         &mut self,
         other_fluid_entity: &mut Self){

        FluidEntityThermophysicalData::
            step_2_conenct_to_component_outlet(
                &mut self.fluid_parameters, 
                &mut other_fluid_entity.fluid_parameters);

    }

    fn step_3_add_component_to_vector(
        &mut self,
        fluid_entity_vector: &mut Vec<FluidEntityThermophysicalData>){

        FluidEntityThermophysicalData::
            step_3_add_component_to_vector(
                &mut self.fluid_parameters, 
                fluid_entity_vector);
    }
}

impl v2_ExplicitCalculationSteps for v2_IterativeHeatFluxTherminolPipe {
    fn step_1_calculate_current_timestep_temp_enthalpies(
        &mut self) {
        // first let's get the temperatures
        // for the "old" inlet and outlet,
        // by old here i mean current timestep
        let current_inlet_temp = 
            self.fluid_parameters.temperature_data.
            inlet_temp_old.clone();

        let current_outlet_temp = 
            self.fluid_parameters.temperature_data.
            outlet_temp_old.clone();

        // now i want to find the average temperature
        // of these temperatures
        // ie T_fluid = (T_in + T_out)/2
        let current_bulk_temp_value_kelvin = 
            (current_inlet_temp.value 
             + current_outlet_temp.value)/2.0;

        let current_bulk_temp = 
            ThermodynamicTemperature::new::<
            thermodynamic_temperature::kelvin>
            (current_bulk_temp_value_kelvin);

        self.fluid_parameters.temperature_data.
            fluid_temp_old = 
            current_bulk_temp.clone();

        // with these values let's find the enthalpies
        self.fluid_parameters.enthalpy_data.
            inlet_enthalpy_old = 
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_inlet_temp);

        self.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old = 
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_outlet_temp);

        self.fluid_parameters.enthalpy_data.
            fluid_enthalpy_old =
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_bulk_temp);
    }

    /// we are calculating
    ///
    /// h_new = h_old + delta T * 
    /// (H_in - H_out + Q + W)
    ///
    /// the fluid properties except for enthalpy
    /// are calculated at T_fluid = (T_in - T_out)/2
    ///
    ///
    /// fluid enthalpy wise, we assume perfect mixing
    /// which means the bulk temperature of the pipe is
    /// the same as the exit temperature
    ///
    /// I haven't quite calculated bidirectional flow yet
    fn step_2_calculate_new_outlet_enthalpy(
        &mut self, 
        heat_supplied_to_fluid: Power,
        work_done_on_fluid: Power,
        timestep: Time,
        fluid_mass_flowrate: MassRate){

        // first let's calculate H_in

        let enthalpy_flowrate_in = 
            fluid_mass_flowrate *
            self.fluid_parameters.
            enthalpy_data.inlet_enthalpy_old.clone();
            
        // second let's calculate H_out

        let enthalpy_flowrate_out = 
            fluid_mass_flowrate *
            self.fluid_parameters.
            enthalpy_data.outlet_enthalpy_old.clone();

        // third let's calculate heat addition

        let heat_addition : Energy = 
            timestep *
            (enthalpy_flowrate_in - enthalpy_flowrate_out 
             + heat_supplied_to_fluid 
             + work_done_on_fluid);

        // fourth, let's calculate the new system enthalpy
        // for this, we need the current of the control
        // volume times specific enthalpy

        let control_volume_mass : Mass = 
            v1_IterativeHeatFluxTherminolPipe::density(
                self.fluid_parameters.temperature_data.
                fluid_temp_old)*
            self.fluid_parameters.fluid_volume;

        let h_old : Energy = 
            control_volume_mass *
            self.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old;

        let h_new : Energy = 
            h_old + heat_addition;

        // fifth, assuming mass does not change much
        // we can find the new outlet specifc enthalpy

        let new_fluid_specifc_enthalpy : AvailableEnergy 
            = h_new/control_volume_mass;

        // now we set the new outlet enthalpy,
        // note that if we want to do a bidirectional flow,
        // we need to do an if else loop here
        // that means, 3 cases
        // what to do if mass flow is +ve,
        // then -ve
        // then 0
        self.fluid_parameters.enthalpy_data.
            outlet_enthalpy_new = new_fluid_specifc_enthalpy;

        // and that's it!
    }

    /// This step calculates the new system temperature
    /// from system specific enthalpy estimates
    fn step_3_calculate_new_outlet_temperature(
        &mut self) -> ThermodynamicTemperature {

        let new_outlet_fluid_temp : ThermodynamicTemperature
            = v1_IterativeHeatFluxTherminolPipe::
            get_temperature_from_enthalpy(
                self.fluid_parameters.enthalpy_data.
                outlet_enthalpy_new);

        self.fluid_parameters.temperature_data.
            outlet_temp_new = 
            new_outlet_fluid_temp.clone();

        return new_outlet_fluid_temp;
    }

    /// This step assumes that the new inlet temperature
    /// has been calculated, and sets the new inlet
    /// temperature as such
    fn step_4_set_inlet_temperature(
        &mut self,
        new_inlet_temperature : ThermodynamicTemperature) {

        self.fluid_parameters.temperature_data.
            inlet_temp_new = new_inlet_temperature;

    }

    /// This step assumes that the new outlet temperature
    /// has been calculated, and sets the new outlet
    /// temperature as such
    fn step_5_set_outlet_temperature(
        &mut self,
        new_outlet_temperature : ThermodynamicTemperature) {

        self.fluid_parameters.temperature_data.
            outlet_temp_new = new_outlet_temperature;

    }


    /// Assuming all new temperatures have been calculated
    /// one can update set the old temperatures
    /// to the new ones at the next timestep
    fn step_6_update_current_timestep_temperatures(
        &mut self) {

        // set inlet and outlet temperatures
        self.fluid_parameters.temperature_data.
            inlet_temp_old =
            self.fluid_parameters.temperature_data.
            inlet_temp_new.clone();

        self.fluid_parameters.temperature_data.
            fluid_temp_old =
            self.fluid_parameters.temperature_data.
            fluid_temp_new.clone();


        // set new fluid temperature to average
        // of new inlet and outlet temperatures

        let new_fluid_temp_value =
            (self.fluid_parameters.temperature_data.
            inlet_temp_new.value + 
            self.fluid_parameters.temperature_data.
            outlet_temp_new.value)/2.0;

        self.fluid_parameters.temperature_data.
            fluid_temp_new = 
            ThermodynamicTemperature::new::
            <kelvin>(new_fluid_temp_value);


        self.fluid_parameters.temperature_data.
            outlet_temp_old =
            self.fluid_parameters.temperature_data.
            outlet_temp_new.clone();

    }
}
use crate::ControlVolumeCalculations::TherminolDowthermPipes::
TherminolFluidProperties;
impl TherminolFluidProperties for v2_IterativeHeatFluxTherminolPipe {
}

/// This struct or class represents a fixed heat
/// flux therminol pipe, or at least the first version
///
///
/// 
/// Here, we don't consider conjugate heat transfer 
/// or anything, but we just supply a fixed heat value to
/// the therminol pipe
///
/// This is the primitive version of calculating 
/// enthalpies and temperatures..
///
/// ```rust
/// 
/// extern crate approx;
/// use heat_transfer_rust::ControlVolumeCalculations::
/// Sandbox::*;
///
/// use uom::si::f64::*;
/// use uom::si::time::second;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use uom::si::volume::cubic_meter;
///
/// // first let's initialise a pipe using a timestep
/// // global thermodynamic temperature and fluid volume
/// 
/// let timestep = Time::new::<second>(0.1_f64);
/// let initial_global_temp = ThermodynamicTemperature::
/// new::<kelvin>(300_f64);
///
/// let fluid_volume = Volume::new::<cubic_meter>(
/// 0.01_f64.powf(3_f64));
///
/// let mut fluid_entity_index: usize = 1;
///
/// // we are now going to initialise stuff
///
/// use heat_transfer_rust::ControlVolumeCalculations::
/// FluidEntity_StructsAndTraits::FluidEntityInitialisationSteps;
///
/// let mut pipe1 = v1_IterativeHeatFluxTherminolPipe::new();
///
/// pipe1.step_0_set_timestep_and_initial_temperatures(
/// timestep,
/// initial_global_temp,
/// fluid_volume,
/// fluid_entity_index);
///
/// // now suppose there are 3 pipes and i want to connect them
///
/// let mut pipe2 = v1_IterativeHeatFluxTherminolPipe::new();
///
/// pipe2.step_0_set_timestep_and_initial_temperatures(
/// timestep,
/// initial_global_temp,
/// fluid_volume,
/// 2);
///
/// let mut pipe3 = v1_IterativeHeatFluxTherminolPipe::new();
///
/// pipe3.step_0_set_timestep_and_initial_temperatures(
/// timestep,
/// initial_global_temp,
/// fluid_volume,
/// 3);
///
///
/// // This section tests if the pipe connection indexing is working 
/// // correctly
/// 
/// // let pipe 2 be connected to the outlet of pipe 1
///
/// pipe1.step_2_conenct_to_component_outlet(&mut pipe2);
///
///
/// // we shall also connect pipe3 to the outlet of pipe2
///
/// pipe2.step_2_conenct_to_component_outlet(&mut pipe3);
///
/// // lastly we shall connect pipe3 to the inlet of pipe1
///
/// pipe1.step_1_connect_to_component_inlet(&mut pipe3);
///
/// // so pipe 1's index should be 1
///
/// assert_eq!(1, pipe1.fluid_parameters.index_data.fluid_entity_index);
///
///
/// // hence the index of the pipe connected to pipe 1's inlet should
/// // be the index of pipe 3 (which is 3)
///
/// assert_eq!(3, pipe1.fluid_parameters.index_data.inlet_fluid_entity_index);
///
/// // finally, the outlet index of pipe 1 should be 2
///
/// assert_eq!(2, pipe1.fluid_parameters.index_data.outlet_fluid_entity_index);
///
///
/// // and for pipe 3, it's index should be 3,
/// // its inlet pipe index should be 2,
/// // and its outlet pipe index should be 1
///
/// assert_eq!(3, pipe3.fluid_parameters.index_data.fluid_entity_index);
/// assert_eq!(2, pipe3.fluid_parameters.index_data.inlet_fluid_entity_index);
/// assert_eq!(1, pipe3.fluid_parameters.index_data.outlet_fluid_entity_index);
///
/// // for pipe 2, the index is 2, pipe 1 is connected to its inlet
/// // and pipe 3 to its outlet
///
/// assert_eq!(2, pipe2.fluid_parameters.index_data.fluid_entity_index);
/// assert_eq!(1, pipe2.fluid_parameters.index_data.inlet_fluid_entity_index);
/// assert_eq!(3, pipe2.fluid_parameters.index_data.outlet_fluid_entity_index);
///
/// // The last part of the tests stipulates that we should 
/// // populate a vector. This vector should be used in a for loop
/// // so that we can autocalculate the temperatures at each time step
/// // 
/// // but for us to know how to calculate things, we'll probably want
/// // to do it manually first, to best know how to operate
/// // I am supposing that pipe 1 has a supply of 100 W
/// // pipe 2 has a heat loss of 20 W and pipe 3 has a heat loss of 80 W
/// // cannot do too much power because it will cause the temperature
/// // to go below 20C, which is out of range
///
/// use uom::si::power::watt;
///
/// let pipe_1_heat = Power::new::<watt>(100_f64);
/// let pipe_2_heat = Power::new::<watt>(-20_f64);
/// let pipe_3_heat = Power::new::<watt>(-80_f64);
///
/// let work_done_on_pipe_rate = Power::new::<watt>(0_f64);
/// 
/// // let's begin step 1 of calculation procedure
///
/// use crate::heat_transfer_rust::ControlVolumeCalculations::
/// ExplicitCalculations::v1_ExplicitCalculationSteps;
///
///
/// pipe1.step_1_calculate_current_timestep_temp_enthalpies();
/// pipe2.step_1_calculate_current_timestep_temp_enthalpies();
/// pipe3.step_1_calculate_current_timestep_temp_enthalpies();
///
/// // then step 2, which is to calculate new system enthalpy
/// // we assume mass flowrate has already been calculated 
/// // for this timestep or is constant
/// 
/// use uom::si::mass_rate::kilogram_per_second;
///
/// let timestep : Time = pipe1.fluid_parameters.timestep;
/// let mass_flowrate : MassRate = MassRate::new::
/// <kilogram_per_second>(0.18);
///
/// pipe1.step_2_calculate_new_system_enthalpy(
/// pipe_1_heat,
/// work_done_on_pipe_rate,
/// timestep,
/// mass_flowrate);
///
///
/// pipe2.step_2_calculate_new_system_enthalpy(
/// pipe_2_heat,
/// work_done_on_pipe_rate,
/// timestep,
/// mass_flowrate);
///
/// 
/// pipe3.step_2_calculate_new_system_enthalpy(
/// pipe_3_heat,
/// work_done_on_pipe_rate,
/// timestep,
/// mass_flowrate);
///
///
///
/// // now for step 3, to calculate new thermodynamic temperature
/// // for now i'll have to debug the therminol properties, it's not
/// // giving the correct number
///
/// let T_new1 = pipe1.step_3_calculate_new_system_temperature();
/// let T_new2 = pipe2.step_3_calculate_new_system_temperature();
/// let T_new3 = pipe3.step_3_calculate_new_system_temperature();
///
/// let mut error_string: String = T_new1.value.to_string();
///
/// error_string.push_str(" ");
/// error_string.push_str(&T_new2.value.to_string());
/// error_string.push_str(" ");
/// error_string.push_str(&T_new3.value.to_string());
///
/// // panic!("{}", error_string);
///
/// // at time of test
/// // T_new1 = 305.91 K
/// // T_new2 = 298.8 K
/// // T_new3 = 295.2 K
///
/// approx::assert_relative_eq!(305.9, T_new1.value, max_relative=0.001);
/// approx::assert_relative_eq!(298.8, T_new2.value, max_relative=0.001);
/// approx::assert_relative_eq!(295.2, T_new3.value, max_relative=0.001);
///
/// //panic!("{}",pipe2.fluid_parameters.temperature_data.
/// //fluid_temp_new.value);
///
/// // to calculate for step 4, i will need to solve a few eqns
/// //  T_new1 (kelvin) = (T_in1 + T_out1)/2
/// //  T_new2 (kelvin) = (T_in2 + T_out2)/2
/// //  T_new3 (kelvin) = (T_in3 + T_out3)/2
/// // 
/// // However, T_in1 = T_out3
/// // T_in2 = T_out_1
/// // T_in3 = T_out2
///
///
/// //  In terms of inlet temperatures, this becomes:
/// //  T_new1 (kelvin) = (T_in1 + T_in2)/2
/// //  T_new2 (kelvin) = (T_in2 + T_in3)/2
/// //  T_new3 (kelvin) = (T_in3 + T_in1)/2
///
/// // if i were to represent this in matrix form
/// //  |0.5 0.5 0.0| | T_in1 |   | T_new1 |
/// //  |0.0 0.5 0.5| | T_in2 | = | T_new2 |
/// //  |0.5 0.0 0.5| | T_in3 |   | T_new3 |
/// //
/// // we then solve this matrix
/// // There are two problems here:
/// // (1) How do we construct this matrix automatically,
/// // (2) How do we solve this matrix automatically
/// // 
///
/// // For now let's solve the second problem first
/// // for reference
/// // solving this results in (i used scilab linsolve)
/// // T_in1 = 302   K
/// // T_in2 = 309.5 K
/// // T_in3 = 288.1 K
/// // 
/// // the way forward:
///
/// extern crate ndarray;
/// use ndarray::prelude::*;
/// // use ndarray_linalg::Solve;
///
/// // first we make a 2D array 
/// let matrixA : Array2<f64> = 
/// array![
/// [0.5, 0.5, 0.0],
/// [0.0, 0.5, 0.5],
/// [0.5, 0.0, 0.5]
/// ];
///
/// let vectorB: Array1<f64> = 
/// array![ 305.91_f64, 298.8, 295.2 ];
/// //
/// // note that for this to work, i needed to install
/// // ndarray_linalg with features
/// // openblas-static
/// // https://github.com/rust-ndarray/ndarray-linalg/issues/171
/// // This means i needed gcc-fortran on my archlinux as well
/// // 
/// // now, openblas is linux only, and intel math kernel
/// // library (intel-mkl) is cross platform but intel only
/// // i am opting not to use openblas to ensure this is cross
/// // platform
///
/// //let x = matrixA.solve_into(vectorB).unwrap();
///
/// //approx::assert_relative_eq!(302_f64, x[0], max_relative=0.01);
/// //approx::assert_relative_eq!(309.5_f64, x[1], max_relative=0.01);
/// //approx::assert_relative_eq!(288.1_f64, x[2], max_relative=0.01);
///
///
/// // now this solves the thing ok, and it means that there is a 7.5 C
/// // increase in pipe 1,
/// // a 21.4C decrease in pipe 2
/// // and a 13.9 C increase in pipe 3
/// // this makes no physical sense
///
/// // as compared to pipe 1 outlet being 305.91K, inlet is 295.2K
/// // pipe 2 outlet 298.8K, inlet 305.91K,
/// // pipe 3 outlet 295.2K, inlet 298.8K 
/// // in the former case, energy is not conserved
/// // therefore, we shall not use the fluid temp
/// // T_fluid = (T_in + T_out)/2.0 
/// // for use in energy balance
/// // thus, i will use the outlet temperature to be the same
/// // temperature as the control volume temperature
/// // similar to a CSTR (continuously stirred tank reactor)
/// // and we assume perfect mixing
/// 
/// // step 4 and 5
/// // for the inlet temperatures we can use the T_new1, T_new2, and T_new3
/// // as outlet temps for pipe 1,2 and 3
///
/// pipe1.step_4_set_inlet_temperature(T_new3);
/// pipe2.step_4_set_inlet_temperature(T_new1);
/// pipe3.step_4_set_inlet_temperature(T_new2);
///
/// pipe1.step_5_set_outlet_temperature(T_new1);
/// pipe2.step_5_set_outlet_temperature(T_new2);
/// pipe3.step_5_set_outlet_temperature(T_new3);
///
/// // step 6: now we set the old inlet and outlet temperatures
///
/// pipe1.step_6_update_current_timestep_temperatures();
/// pipe2.step_6_update_current_timestep_temperatures();
/// pipe3.step_6_update_current_timestep_temperatures();
///
/// // if everything works out well, the old temperatures
/// // should reflect the new values at the next time step
///
///
/// approx::assert_relative_eq!(305.9, 
/// pipe1.fluid_parameters.temperature_data.outlet_temp_old.value, 
/// max_relative=0.001);
///
/// approx::assert_relative_eq!(298.8, 
/// pipe2.fluid_parameters.temperature_data.outlet_temp_old.value, 
/// max_relative=0.001);
///
/// approx::assert_relative_eq!(295.2, 
/// pipe3.fluid_parameters.temperature_data.outlet_temp_old.value, 
/// max_relative=0.001);
///
/// // now we test the inlet temperatures
///
/// approx::assert_relative_eq!(305.9, 
/// pipe2.fluid_parameters.temperature_data.inlet_temp_old.value, 
/// max_relative=0.001);
///
/// approx::assert_relative_eq!(298.8, 
/// pipe3.fluid_parameters.temperature_data.inlet_temp_old.value, 
/// max_relative=0.001);
///
/// approx::assert_relative_eq!(295.2, 
/// pipe1.fluid_parameters.temperature_data.inlet_temp_old.value, 
/// max_relative=0.001);
/// 
/// // this concludes the testing for one round of a 3 pipe problem
/// // we might want to do the temperature assignment in a more streamlined
/// // manner
///
/// ```
///
///
///
///
///
#[derive(Clone)]
pub struct v1_IterativeHeatFluxTherminolPipe {
    pub fluid_parameters: FluidEntityThermophysicalData,
}

impl v1_IterativeHeatFluxTherminolPipe {

    /// Constructor which creates the structure 
    /// ```rust
    /// 
    /// extern crate approx;
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// TherminolDowthermPipes::*;
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// Sandbox::*;
    ///
    /// use uom::si::f64::*;
    /// use uom::si::time::second;
    /// use uom::si::thermodynamic_temperature::kelvin;
    /// use uom::si::volume::cubic_meter;
    ///
    /// // first let's initialise a pipe using a timestep
    /// // global thermodynamic temperature and fluid volume
    /// 
    /// let timestep = Time::new::<second>(0.1_f64);
    /// let initial_global_temp = ThermodynamicTemperature::
    /// new::<kelvin>(300_f64);
    ///
    /// let fluid_volume = Volume::new::<cubic_meter>(
    /// 0.01_f64.powf(3_f64));
    ///
    /// let fluid_entity_index: usize = 4;
    ///
    /// // we are now going to initialise stuff
    ///
    /// use heat_transfer_rust::ControlVolumeCalculations::
    /// FluidEntity_StructsAndTraits::FluidEntityInitialisationSteps;
    ///
    /// let mut pipe1 = v1_IterativeHeatFluxTherminolPipe::new();
    ///
    /// pipe1.step_0_set_timestep_and_initial_temperatures(
    /// timestep,
    /// initial_global_temp,
    /// fluid_volume,
    /// fluid_entity_index);
    ///
    ///
    ///
    /// ```
    pub fn new() -> Self {

        // we'll need to make a few data structures first
        // with some default values


        let default_timestep: Time = 
            Time::new::<second>(0.1);

        let default_fluid_volume: Volume = 
            Volume::new::<cubic_meter>(0.02_f64.powf(3.0));

        // let's populate default index data
        
        let default_index: usize = 0;

        let default_index_data = 
            FluidEntityIndexData { 
                fluid_entity_index : default_index ,
                inlet_fluid_entity_index : default_index,
                outlet_fluid_entity_index : default_index,
            };


        // let's populate default temperature data
        //

        let default_temperature : ThermodynamicTemperature 
            = ThermodynamicTemperature::new::
            <kelvin>(310_f64);

        let default_temperature_data : PipeFluidTemperatureData  
            = PipeFluidTemperatureData { 
                inlet_temp_old: default_temperature, 
                inlet_temp_new: default_temperature, 
                outlet_temp_old: default_temperature, 
                outlet_temp_new: default_temperature, 
                fluid_temp_old: default_temperature, 
                fluid_temp_new: default_temperature,
            };

        // let's finally populate the enthalpy data

        let default_enthalpy : AvailableEnergy =
            AvailableEnergy::new::<joule_per_kilogram>(0_f64);

        let default_enthalpy_data : PipeFluidEnthalpyData =
            PipeFluidEnthalpyData { 
                inlet_enthalpy_old: default_enthalpy, 
                inlet_enthalpy_new: default_enthalpy, 
                outlet_enthalpy_old: default_enthalpy, 
                outlet_enthalpy_new: default_enthalpy, 
                fluid_enthalpy_old: default_enthalpy, 
                fluid_enthalpy_new: default_enthalpy 
            };

        let default_thermophysical_data : FluidEntityThermophysicalData
            = FluidEntityThermophysicalData { 
                index_data: default_index_data, 
                temperature_data: default_temperature_data, 
                enthalpy_data: default_enthalpy_data, 
                timestep: default_timestep, 
                fluid_volume: default_fluid_volume ,
            };

        return Self { 
            fluid_parameters: default_thermophysical_data 
        };

    }
}

impl FluidEntityInitialisationSteps for v1_IterativeHeatFluxTherminolPipe {
    /// step zero, essentially the constructor

    fn step_0_set_timestep_and_initial_temperatures(
        &mut self,
        timestep: Time,
        initial_global_temp: ThermodynamicTemperature,
        fluid_volume: Volume,
        fluid_entity_index: usize) -> Self {

        FluidEntityThermophysicalData::
            step_0_set_timestep_and_initial_temperatures(
                &mut self.fluid_parameters, 
                timestep, 
                initial_global_temp, 
                fluid_volume, 
                fluid_entity_index
                );

        return self.clone();
    }
    

    fn step_1_connect_to_component_inlet(
        &mut self,
        other_fluid_entity: &mut Self) {

        FluidEntityThermophysicalData::
            step_1_connect_to_component_inlet(
                &mut self.fluid_parameters, 
                &mut other_fluid_entity.fluid_parameters);
    }

    fn step_2_conenct_to_component_outlet(
         &mut self,
         other_fluid_entity: &mut Self){

        FluidEntityThermophysicalData::
            step_2_conenct_to_component_outlet(
                &mut self.fluid_parameters, 
                &mut other_fluid_entity.fluid_parameters);

    }

    fn step_3_add_component_to_vector(
        &mut self,
        fluid_entity_vector: &mut Vec<FluidEntityThermophysicalData>){

        FluidEntityThermophysicalData::
            step_3_add_component_to_vector(
                &mut self.fluid_parameters, 
                fluid_entity_vector);
    }
}

impl v1_ExplicitCalculationSteps for v1_IterativeHeatFluxTherminolPipe {
    fn step_1_calculate_current_timestep_temp_enthalpies(
        &mut self) {
        // first let's get the temperatures
        //
        let current_inlet_temp = 
            self.fluid_parameters.temperature_data.
            inlet_temp_old.clone();

        let current_outlet_temp = 
            self.fluid_parameters.temperature_data.
            outlet_temp_old.clone();

        // now i want to find the average temperature
        // of these temperatures
        // ie T_fluid = (T_in + T_out)/2
        let current_bulk_temp_value_kelvin = 
            (current_inlet_temp.value 
             + current_outlet_temp.value)/2.0;

        let current_bulk_temp = 
            ThermodynamicTemperature::new::<
            thermodynamic_temperature::kelvin>
            (current_bulk_temp_value_kelvin);

        self.fluid_parameters.temperature_data.
            fluid_temp_old = 
            current_bulk_temp.clone();

        // with these values let's find the enthalpies
        self.fluid_parameters.enthalpy_data.
            inlet_enthalpy_old = 
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_inlet_temp);

        self.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old = 
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_outlet_temp);

        self.fluid_parameters.enthalpy_data.
            fluid_enthalpy_old =
            v1_IterativeHeatFluxTherminolPipe::
            enthalpy(current_bulk_temp);
    }

    /// we are calculating
    ///
    /// h_new = h_old + delta T * 
    /// (H_in - H_out + Q + W)
    ///
    /// the fluid properties except for enthalpy
    /// are calculated at T_fluid = (T_in - T_out)/2
    ///
    ///
    /// fluid enthalpy wise, we assume perfect mixing
    /// which means the bulk temperature of the pipe is
    /// the same as the exit temperature
    fn step_2_calculate_new_system_enthalpy(
        &mut self, 
        heat_supplied_to_fluid: Power,
        work_done_on_fluid: Power,
        timestep: Time,
        fluid_mass_flowrate: MassRate){

        // first let's calculate H_in

        let enthalpy_flowrate_in = 
            fluid_mass_flowrate *
            self.fluid_parameters.
            enthalpy_data.inlet_enthalpy_old.clone();
            
        // second let's calculate H_out

        let enthalpy_flowrate_out = 
            fluid_mass_flowrate *
            self.fluid_parameters.
            enthalpy_data.outlet_enthalpy_old.clone();

        // third let's calculate heat addition

        let heat_addition : Energy = 
            timestep *
            (enthalpy_flowrate_in - enthalpy_flowrate_out 
             + heat_supplied_to_fluid 
             + work_done_on_fluid);

        // fourth, let's calculate the new system enthalpy
        // for this, we need the current of the control
        // volume times specific enthalpy

        let control_volume_mass : Mass = 
            v1_IterativeHeatFluxTherminolPipe::density(
                self.fluid_parameters.temperature_data.
                fluid_temp_old)*
            self.fluid_parameters.fluid_volume;

        let h_old : Energy = 
            control_volume_mass *
            self.fluid_parameters.enthalpy_data.
            outlet_enthalpy_old;

        let h_new : Energy = 
            h_old + heat_addition;

        // fifth, assuming mass does not change much
        // we can find the new system specifc enthalpy

        let new_fluid_specifc_enthalpy : AvailableEnergy 
            = h_new/control_volume_mass;

        self.fluid_parameters.enthalpy_data.
            fluid_enthalpy_new = new_fluid_specifc_enthalpy;

        // and that's it!
    }

    /// This step calculates the new system temperature
    /// from system specific enthalpy estimates
    fn step_3_calculate_new_system_temperature(
        &mut self) -> ThermodynamicTemperature {

        let new_fluid_temp_average : ThermodynamicTemperature
            = v1_IterativeHeatFluxTherminolPipe::
            get_temperature_from_enthalpy(
                self.fluid_parameters.enthalpy_data.
                fluid_enthalpy_new);

        self.fluid_parameters.temperature_data.
            fluid_temp_new = 
            new_fluid_temp_average.clone();

        return new_fluid_temp_average;
    }

    /// This step assumes that the new inlet temperature
    /// has been calculated, and sets the new inlet
    /// temperature as such
    fn step_4_set_inlet_temperature(
        &mut self,
        new_inlet_temperature : ThermodynamicTemperature) {

        self.fluid_parameters.temperature_data.
            inlet_temp_new = new_inlet_temperature;

    }

    /// This step assumes that the new outlet temperature
    /// has been calculated, and sets the new outlet
    /// temperature as such
    fn step_5_set_outlet_temperature(
        &mut self,
        new_outlet_temperature : ThermodynamicTemperature) {

        self.fluid_parameters.temperature_data.
            outlet_temp_new = new_outlet_temperature;

    }


    /// Assuming all new temperatures have been calculated
    /// one can update set the old temperatures
    /// to the new ones at the next timestep
    fn step_6_update_current_timestep_temperatures(
        &mut self) {

        // set inlet and outlet temperatures
        self.fluid_parameters.temperature_data.
            inlet_temp_old =
            self.fluid_parameters.temperature_data.
            inlet_temp_new.clone();

        self.fluid_parameters.temperature_data.
            fluid_temp_old =
            self.fluid_parameters.temperature_data.
            fluid_temp_new.clone();


        // set new fluid temperature to average
        // of new inlet and outlet temperatures

        let new_fluid_temp_value =
            (self.fluid_parameters.temperature_data.
            inlet_temp_new.value + 
            self.fluid_parameters.temperature_data.
            outlet_temp_new.value)/2.0;

        self.fluid_parameters.temperature_data.
            fluid_temp_new = 
            ThermodynamicTemperature::new::
            <kelvin>(new_fluid_temp_value);


        self.fluid_parameters.temperature_data.
            outlet_temp_old =
            self.fluid_parameters.temperature_data.
            outlet_temp_new.clone();

    }
}
impl TherminolFluidProperties for v1_IterativeHeatFluxTherminolPipe {
}
