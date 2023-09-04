
use crate::{fluid_mechanics_lib::fluid_component_calculation::fluid_component_trait, heat_transfer_lib::thermophysical_properties::{dynamic_viscosity::dynamic_viscosity, density::density}};
use fluid_component_trait::FluidComponent;

use super::FluidArray;
use uom::si::f64::*;

impl FluidComponent for FluidArray {
    fn get_mass_flowrate(&mut self) -> MassRate  {
        self.mass_flowrate
    }

    fn set_mass_flowrate(&mut self, mass_flowrate: MassRate) {
        self.mass_flowrate = mass_flowrate
    }

    fn get_mass_flowrate_from_pressure_loss_immutable(
        &self, pressure_loss: Pressure) -> MassRate {
        todo!("need pressure loss to Re Correlatoin")
    }

    fn get_pressure_loss(&mut self) -> Pressure {
        self.pressure_loss
    }

    fn set_pressure_loss(&mut self, pressure_loss: Pressure) {
        self.pressure_loss = pressure_loss
    }

    fn get_pressure_loss_immutable(
        &self, mass_flowrate: MassRate) -> Pressure {
        todo!("need pressure loss to Re Correlatoin")
    }

    fn get_cross_sectional_area(&mut self) -> Area {
        self.xs_area
    }

    fn get_cross_sectional_area_immutable(&self) -> Area {
        self.xs_area
    }

    fn get_hydraulic_diameter(&mut self) -> Length {
        // d_h = 4A/P 

        4.0 * self.xs_area / self.wetted_perimiter
    }

    fn get_hydraulic_diameter_immutable(&self) -> Length {
        4.0 * self.xs_area / self.wetted_perimiter
    }

    fn get_fluid_viscosity(&mut self) -> DynamicViscosity {
        let temperature = self.get_bulk_temperature().unwrap();

        let viscosity = dynamic_viscosity(
            self.material_control_volume,
            temperature,
            self.pressure_control_volume).unwrap();

        return viscosity;
    }

    fn get_fluid_viscosity_immutable(&self) -> DynamicViscosity {
        let temperature = self.clone().get_bulk_temperature().unwrap();

        let viscosity = dynamic_viscosity(
            self.material_control_volume,
            temperature,
            self.pressure_control_volume).unwrap();

        return viscosity;
    }

    fn get_fluid_density(&mut self) -> MassDensity {
        let temperature = self.get_bulk_temperature().unwrap();

        let density = density(
            self.material_control_volume,
            temperature,
            self.pressure_control_volume).unwrap();

        return density;
    }

    fn get_fluid_density_immutable(&self) -> MassDensity {
        let temperature = self.clone().get_bulk_temperature().unwrap();

        let density = density(
            self.material_control_volume,
            temperature,
            self.pressure_control_volume).unwrap();

        return density;
    }

    fn get_component_length(&mut self) -> Length {
        self.total_length
    }

    fn get_component_length_immutable(&self) -> Length {
        self.total_length
    }

    fn get_incline_angle(&mut self) -> Angle {
        self.incline_angle
    }

    fn get_incline_angle_immutable(&self) -> Angle {
        self.incline_angle
    }

    fn get_internal_pressure_source(&mut self) -> Pressure {
        self.internal_pressure_source
    }

    fn get_internal_pressure_source_immutable(&self) -> Pressure {
        self.internal_pressure_source
    }

    fn set_internal_pressure_source(
        &mut self,
        internal_pressure: Pressure) {
        self.internal_pressure_source = internal_pressure;
    }
}
