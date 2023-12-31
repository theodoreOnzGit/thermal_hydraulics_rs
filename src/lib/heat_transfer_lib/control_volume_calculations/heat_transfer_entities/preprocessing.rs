use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::link_heat_transfer_entity;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::HeatTransferInteractionType;
use crate::fluid_mechanics_lib::prelude::FluidComponent;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::heat_transfer_lib::thermophysical_properties::Material;

use super::{HeatTransferEntity, CVType};
use uom::si::f64::*;

impl HeatTransferEntity {

    /// get maximum timestep 
    ///
    /// 
    pub fn get_max_timestep(
        entity: &mut HeatTransferEntity,
        max_temperature_change: TemperatureInterval) 
    -> Result<Time, String> {

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => 
                return Err("getting timestep not \n 
                    implemented for BoundaryConditions".to_string()),
        };

        let cv_timestep_result: 
        Result<Time, String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.get_max_timestep(max_temperature_change)
                },
                CVType::ArrayCV(cv) => {
                    cv.get_max_timestep(max_temperature_change)
                },
            };

        return cv_timestep_result;
    }

    /// tries to set or get mass flowrate from the heat transfer 
    /// entity 
    ///
    /// only works in the case of the fluid array 
    /// everything else does not matter
    pub fn try_set_mass_flowrate(
    entity: &mut HeatTransferEntity,
    mass_flowrate: MassRate) -> Result<(),ThermalHydraulicsLibError> {

        match entity {
            HeatTransferEntity::ControlVolume(cv_type) => {
                match cv_type {
                    CVType::ArrayCV(array_cv) => {
                        match array_cv {
                            super::ArrayCVType::GenericPipe(fluid_arr) => {
                                fluid_arr.set_mass_flowrate(mass_flowrate)
                            },
                            _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
                            heat transfer entity"),
                        }
                    },
                    _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
                    heat transfer entity"),
                }
            },
            _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
            heat transfer entity"),
        }

        Ok(())
    }
    /// tries to set or get mass flowrate from the heat transfer 
    /// entity 
    ///
    /// only works in the case of the fluid array 
    /// everything else does not matter
    pub fn get_mass_flowrate(
    entity: &mut HeatTransferEntity) 
    -> Result<MassRate,ThermalHydraulicsLibError> {

        match entity {
            HeatTransferEntity::ControlVolume(cv_type) => {
                match cv_type {
                    CVType::ArrayCV(array_cv) => {
                        match array_cv {
                            super::ArrayCVType::GenericPipe(fluid_arr) => {
                                Ok(fluid_arr.get_mass_flowrate())
                            },
                            _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
                            heat transfer entity"),
                        }
                    },
                    _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
                    heat transfer entity"),
                }
            },
            _ => unimplemented!("mass flowrate irrelevant for this kind of  \n
            heat transfer entity"),
        }

    }

    /// wrapper for linking, makes it easier to link 
    #[inline]
    pub fn link_to_front(&mut self,
    other_hte: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError>{

        link_heat_transfer_entity(
            self,
            other_hte,
            interaction).unwrap();

        Ok(())
    }

    /// wrapper for linking, makes it easier to link 
    #[inline]
    pub fn link_to_back(&mut self,
    other_hte: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError>{

        link_heat_transfer_entity(
            other_hte,
            self,
            interaction).unwrap();

        Ok(())
    }

    /// wrapper for linking, makes it easier to link 
    #[inline]
    pub fn link(entity: &mut HeatTransferEntity,
    other_hte: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError>{

        link_heat_transfer_entity(
            entity,
            other_hte,
            interaction).unwrap();

        Ok(())
    }

    /// tries to get the material of the control volume 
    #[inline]
    pub fn try_get_material_control_volume(&mut self) 
        -> Result<Material, ThermalHydraulicsLibError> {


        let mut cv: CVType = self.clone().try_into()?;

        cv.get_material()

    }
}
