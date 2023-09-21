use thermal_hydraulics_rs::prelude::alpha_nightly::*;
use super::HeaterVersion2Bare;


impl FluidComponent for HeaterVersion2Bare {
    fn get_mass_flowrate(&mut self) -> MassRate  {
        let mut therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_mass_flowrate()
    }

    fn set_mass_flowrate(&mut self, mass_flowrate: MassRate) {
        
        // the therminol array must be accessed first 
        // and then cloned

        let mut therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.set_mass_flowrate(mass_flowrate);

        // unfortunately, this makes setting mass flowrate quite 
        // expensive as we need to clone it everytime

        self.therminol_array.set(therminol_array.into()).unwrap();
    }

    fn get_mass_flowrate_from_pressure_loss_immutable(
        &self, pressure_loss: Pressure) -> MassRate {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_mass_flowrate_from_pressure_loss_immutable(
            pressure_loss)
    }

    fn get_pressure_loss(&mut self) -> Pressure {
        let mut therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_pressure_loss()
    }

    fn set_pressure_loss(&mut self, pressure_loss: Pressure) {
        let mut therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.set_pressure_loss(pressure_loss)
    }

    fn get_pressure_loss_immutable(
        &self, mass_flowrate: MassRate) -> Pressure {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_pressure_loss_immutable(mass_flowrate)
    }

    fn get_cross_sectional_area(&mut self) -> Area {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_cross_sectional_area_immutable()
    }

    fn get_cross_sectional_area_immutable(&self) -> Area {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_cross_sectional_area_immutable()
    }

    fn get_hydraulic_diameter(&mut self) -> Length {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_hydraulic_diameter_immutable()
    }

    fn get_hydraulic_diameter_immutable(&self) -> Length {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_hydraulic_diameter_immutable()
    }

    fn get_fluid_viscosity(&mut self) -> DynamicViscosity {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_fluid_viscosity_immutable()
    }

    fn get_fluid_viscosity_immutable(&self) -> DynamicViscosity {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_fluid_viscosity_immutable()
    }

    fn get_fluid_density(&mut self) -> MassDensity {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_fluid_density_immutable()
    }

    fn get_fluid_density_immutable(&self) -> MassDensity {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_fluid_density_immutable()
    }

    fn get_component_length(&mut self) -> Length {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_component_length_immutable()
    }

    fn get_component_length_immutable(&self) -> Length {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_component_length_immutable()
    }

    fn get_incline_angle(&mut self) -> Angle {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_incline_angle_immutable()
    }

    fn get_incline_angle_immutable(&self) -> Angle {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_incline_angle_immutable()
    }

    fn get_internal_pressure_source(&mut self) -> Pressure {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_internal_pressure_source_immutable()
    }

    fn get_internal_pressure_source_immutable(&self) -> Pressure {
        let therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.get_internal_pressure_source_immutable()
    }

    fn set_internal_pressure_source(
        &mut self,
        internal_pressure: Pressure) {
        let mut therminol_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        therminol_array.set_pressure_loss(internal_pressure)
    }
}
