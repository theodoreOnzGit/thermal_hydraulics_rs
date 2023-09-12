use thermal_hydraulics_rs::prelude::alpha_nightly::*;
use super::HeaterVersion2Bare;


impl FluidComponent for HeaterVersion2Bare {
    fn get_mass_flowrate(&mut self) -> MassRate  {
        todo!()
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
        todo!()
    }

    fn get_pressure_loss(&mut self) -> Pressure {
        todo!()
    }

    fn set_pressure_loss(&mut self, pressure_loss: Pressure) {
        todo!()
    }

    fn get_pressure_loss_immutable(
        &self, mass_flowrate: MassRate) -> Pressure {
        todo!()
    }

    fn get_cross_sectional_area(&mut self) -> Area {
        todo!()
    }

    fn get_cross_sectional_area_immutable(&self) -> Area {
        todo!()
    }

    fn get_hydraulic_diameter(&mut self) -> Length {
        todo!()
    }

    fn get_hydraulic_diameter_immutable(&self) -> Length {
        todo!()
    }

    fn get_fluid_viscosity(&mut self) -> DynamicViscosity {
        todo!()
    }

    fn get_fluid_viscosity_immutable(&self) -> DynamicViscosity {
        todo!()
    }

    fn get_fluid_density(&mut self) -> MassDensity {
        todo!()
    }

    fn get_fluid_density_immutable(&self) -> MassDensity {
        todo!()
    }

    fn get_component_length(&mut self) -> Length {
        todo!()
    }

    fn get_component_length_immutable(&self) -> Length {
        todo!()
    }

    fn get_incline_angle(&mut self) -> Angle {
        todo!()
    }

    fn get_incline_angle_immutable(&self) -> Angle {
        todo!()
    }

    fn get_internal_pressure_source(&mut self) -> Pressure {
        todo!()
    }

    fn get_internal_pressure_source_immutable(&self) -> Pressure {
        todo!()
    }

    fn set_internal_pressure_source(
        &mut self,
        internal_pressure: Pressure) {
        todo!()
    }
}
