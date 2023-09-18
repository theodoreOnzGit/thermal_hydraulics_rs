use uom::{si::f64::*, num_traits::Zero};

use crate::prelude::alpha_nightly::ThermalHydraulicsLibError;

use super::HeatTransferEntity;
/// Contains all the types of Boundary Conditions (BCs) you can use 
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum BCType {
    /// The user specifies a fixed temperature for the BC
    UserSpecifiedTemperature(ThermodynamicTemperature),
    /// The user specifies a heat flux for the BC
    /// the uom type is heat flux density in power/area
    UserSpecifiedHeatFlux(HeatFluxDensity),
    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}

impl BCType {
    /// creates a new constant temperature BC
    pub fn new_const_temperature(temperature:ThermodynamicTemperature)
        -> HeatTransferEntity {
        return HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedTemperature(temperature));
    }

    /// creates a new constant heat flux bc
    pub fn new_const_heat_flux(heat_flux: HeatFluxDensity)
        -> HeatTransferEntity {
        return HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatFlux(heat_flux));
    }

    /// creates a new constant heat addition bc
    pub fn new_const_heat_addition(heat_addition: Power)
        -> HeatTransferEntity {
        return HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatAddition(heat_addition));
    }

    /// creates a new constant heat addition bc
    pub fn new_adiabatic_bc() -> HeatTransferEntity {
        return HeatTransferEntity::BoundaryConditions(
        BCType::UserSpecifiedHeatAddition(
            Power::zero()));
    }

    pub(crate) fn get_temperature_vector(&self) -> 
    Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{

        match self {
            BCType::UserSpecifiedTemperature(temperature) => {
                let mut temp_vec = vec![];
                temp_vec.push(*temperature);
                return Ok(temp_vec);
            },
            BCType::UserSpecifiedHeatFlux(_) => unimplemented!(),
            BCType::UserSpecifiedHeatAddition(_) => unimplemented!(),
        }
    }


}

/// converts BC types into heat transfer entities and back
pub mod type_conversion;
pub use type_conversion::*;
