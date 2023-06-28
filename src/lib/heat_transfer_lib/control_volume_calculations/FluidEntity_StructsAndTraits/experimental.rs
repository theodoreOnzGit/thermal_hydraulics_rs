extern crate uom;
use uom::si::f64::*;
/// This struct contains information for Fluid
/// Temperatures in a Pipe (ie one inlet and
/// one outlet)
#[derive(Clone, Debug)]
pub struct PipeFluidTemperatureData {
    pub inlet_temp_old: ThermodynamicTemperature,
    pub inlet_temp_new: ThermodynamicTemperature,
    pub outlet_temp_old: ThermodynamicTemperature,
    pub outlet_temp_new: ThermodynamicTemperature,
    pub fluid_temp_old: ThermodynamicTemperature,
    pub fluid_temp_new: ThermodynamicTemperature
}

// to do: focus more on trait objects

/// This struct contains information for
/// Fluid Enthalpy in a Pipe, ie inlet and 
/// outlet
#[derive(Clone)]
pub struct PipeFluidEnthalpyData {
    pub inlet_enthalpy_old: AvailableEnergy,
    pub inlet_enthalpy_new: AvailableEnergy,
    pub outlet_enthalpy_old: AvailableEnergy,
    pub outlet_enthalpy_new: AvailableEnergy,
    pub fluid_enthalpy_old: AvailableEnergy,
    pub fluid_enthalpy_new: AvailableEnergy,

}

/// This structure stores the index
/// of the fluid entity (pipe or some other component)
/// 
/// as well as the indices of the pipes or fluid entities
/// connected to the inlet and outlet
#[derive(Clone)]
pub struct FluidEntityIndexData {
    pub fluid_entity_index: usize,
    pub inlet_fluid_entity_index: usize,
    pub outlet_fluid_entity_index: usize,
}

/// This structure stores the basic data for a 
/// fluid entity
/// 
#[derive(Clone)]
pub struct FluidEntityThermophysicalData {
    pub index_data: FluidEntityIndexData,
    pub temperature_data: PipeFluidTemperatureData,
    pub enthalpy_data: PipeFluidEnthalpyData,
    pub timestep: Time,
    pub fluid_volume: Volume
}




pub trait FluidEntityInitialisationSteps {

    /// Step zero: set timestep and initial temperautres
    ///
    /// Also, the fluid volume for the fluid portion of the
    /// pipe can be assumed fixed (in this case we ignore 
    /// thermal expansion for simplicity)
    /// Otherwise, fluid volume and fluid density must be
    /// taken at each timestep as appropriate parameters
    fn step_0_set_timestep_and_initial_temperatures(
        &mut self,
        timestep: Time,
        initial_global_temp: ThermodynamicTemperature,
        fluid_volume: Volume,
        fluid_entity_index: usize) -> Self;

    /// Step 1: connect a pipe or some other structure
    /// to the inlet to this component or fluid entity
    fn step_1_connect_to_component_inlet(
        &mut self,
        other_fluid_entity: &mut Self);

    /// Step 2: connect a pipe or some other structure
    /// to the outlet of this component or fluid entity
    ///
    /// This step is optional because step 1 should be
    /// able to connect pipe A's inlet to pipe B's outlet
    fn step_2_conenct_to_component_outlet(
        &mut self,
        other_fluid_entity: &mut Self);

    /// Step 3: add component to list or vector of components
    fn step_3_add_component_to_vector(
        &mut self,
        fluid_entity_vector: &mut Vec<FluidEntityThermophysicalData>
        );

}

impl FluidEntityInitialisationSteps 
for FluidEntityThermophysicalData {
    fn step_0_set_timestep_and_initial_temperatures(
        &mut self,
        timestep: Time,
        initial_global_temp: ThermodynamicTemperature,
        fluid_volume: Volume,
        fluid_entity_index: usize) -> Self {
        
        self.fluid_volume = fluid_volume;


        // set temperatures

        self.temperature_data.inlet_temp_old = initial_global_temp;
        self.temperature_data.inlet_temp_new = initial_global_temp;
        self.temperature_data.outlet_temp_old = initial_global_temp;
        self.temperature_data.outlet_temp_new = initial_global_temp;
        self.temperature_data.fluid_temp_old = initial_global_temp;
        self.temperature_data.fluid_temp_new = initial_global_temp;


        // assign index number
        self.index_data.fluid_entity_index = fluid_entity_index;

        // return self
        return self.clone();
    }

    fn step_1_connect_to_component_inlet(
        &mut self,
        other_fluid_entity: &mut FluidEntityThermophysicalData){

        // to connect another component to this component's
        // inlet, i recognise that i must assign the other
        // fluid component's index number to the 
        // index number of the outlet for self

        self.index_data.inlet_fluid_entity_index = 
            other_fluid_entity.index_data.fluid_entity_index.
            clone();
        
        // I likewise ensure that the other fluid entity
        // index has its outlet fluid entity index
        // assigned to this fluid entity index

        other_fluid_entity.index_data.
            outlet_fluid_entity_index = 
            self.index_data.fluid_entity_index.clone();

    }

    fn step_2_conenct_to_component_outlet(
        &mut self,
        other_fluid_entity: &mut FluidEntityThermophysicalData){

        // for this, to connect another component
        // to this component's outlet,  i look for the
        // other fluid_entity's index and then 
        // look for it's outlet index
        
        self.index_data.outlet_fluid_entity_index = 
            other_fluid_entity.index_data.fluid_entity_index.
            clone();

        // i then connect this component to the other 
        // component's inlet
        //
        other_fluid_entity.index_data.
            inlet_fluid_entity_index =
            self.index_data.fluid_entity_index.clone();

    }

    fn step_3_add_component_to_vector(
        &mut self,
        fluid_entity_vector: &mut Vec<FluidEntityThermophysicalData>
        ){

        // here i just push a clone of the fluid entity
        // into the vector
        // it need not be arranged in any particular
        // order
        fluid_entity_vector.push(self.clone());
    }


}


