use crate::boussinesq_solver::array_control_vol::one_d_fluid_array_with_lateral_coupling::FluidArray;
use uom::si::f64::*;

use super::fluid_component_trait::{self, FluidComponentTrait};


/// FluidComponents are pipes and fittings you can connect in parallel
/// such that you can calculate mass flowrate and pressure drop from them
pub enum FluidComponent {
    /// these are arrays of control volumes connected in series
    FluidArray(FluidArray),
}

impl FluidComponentTrait for FluidComponent {
    fn get_mass_flowrate(&mut self) -> MassRate  {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_mass_flowrate()
            },
        }
    }

    fn set_mass_flowrate(&mut self, mass_flowrate: MassRate) {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.set_mass_flowrate(mass_flowrate)
            },
        }
    }

    fn get_mass_flowrate_from_pressure_loss_immutable(
        &self, pressure_loss: Pressure) -> MassRate {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_mass_flowrate_from_pressure_loss_immutable(pressure_loss)
            },
        }
    }

    fn get_pressure_loss(&mut self) -> Pressure {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_pressure_loss()
            },
        }
    }

    fn set_pressure_loss(&mut self, pressure_loss: Pressure) {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.set_pressure_loss(pressure_loss)
            },
        }
    }

    fn get_pressure_loss_immutable(
        &self, mass_flowrate: MassRate) -> Pressure {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_pressure_loss_immutable(mass_flowrate)
            },
        }
    }

    fn get_cross_sectional_area(&mut self) -> Area {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_cross_sectional_area()
            },
        }
    }

    fn get_cross_sectional_area_immutable(&self) -> Area {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_cross_sectional_area_immutable()
            },
        }
    }

    fn get_hydraulic_diameter(&mut self) -> Length {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_hydraulic_diameter()
            },
        }
    }

    fn get_hydraulic_diameter_immutable(&self) -> Length {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_hydraulic_diameter_immutable()
            },
        }
    }

    fn get_fluid_viscosity(&mut self) -> DynamicViscosity {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_fluid_viscosity()
            },
        }
    }

    fn get_fluid_viscosity_immutable(&self) -> DynamicViscosity {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_fluid_viscosity_immutable()
            },
        }
    }

    fn get_fluid_density(&mut self) -> MassDensity {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_fluid_density()
            },
        }
    }

    fn get_fluid_density_immutable(&self) -> MassDensity {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_fluid_density_immutable()
            },
        }
    }

    fn get_component_length(&mut self) -> Length {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_component_length()
            },
        }
    }

    fn get_component_length_immutable(&self) -> Length {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_component_length_immutable()
            },
        }
    }

    fn get_incline_angle(&mut self) -> Angle {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_incline_angle()
            },
        }
    }

    fn get_incline_angle_immutable(&self) -> Angle {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_incline_angle_immutable()
            },
        }
    }

    fn get_internal_pressure_source(&mut self) -> Pressure {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_internal_pressure_source()
            },
        }
    }

    fn get_internal_pressure_source_immutable(&self) -> Pressure {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.get_internal_pressure_source_immutable()
            },
        }
    }

    fn set_internal_pressure_source(
        &mut self,
        internal_pressure: Pressure) {
        match self {
            FluidComponent::FluidArray(fluid_array) => {
                fluid_array.set_internal_pressure_source(internal_pressure)
            },
        }
    }
}
