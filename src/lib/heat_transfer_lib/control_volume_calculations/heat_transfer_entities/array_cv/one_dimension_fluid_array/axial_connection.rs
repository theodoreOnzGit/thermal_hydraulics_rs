use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::HeatTransferEntity;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::HeatTransferInteractionType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::link_heat_transfer_entity;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_entities::CVType;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::FluidArray;

/// this implementation deals with lateral connections 
///
/// the convention is to supply an average conductance 
/// as well as a temperature array
///
/// at the end of the connection phase, one can then use 
/// the advance_timestep method to calculate the new 
/// temperature array
impl FluidArray{

    /// connects an adjacent heat transfer entity to front cv 
    ///
    pub fn link_heat_transfer_entity_to_front_cv(&mut self,
    heat_trf_entity: &mut HeatTransferEntity,
    heat_transfer_interaction: HeatTransferInteractionType) 
        -> Result<(), ThermalHydraulicsLibError>{

        // let's find the front cv, 
        // we shall clone it, link it, then replace the front cv 
        // with the edited control volume 

        let front_cv_clone = self.front_single_cv.clone();

        // package it as a heat transfer entity 

        let mut front_cv: HeatTransferEntity = 
        HeatTransferEntity::ControlVolume(
            CVType::SingleCV(front_cv_clone)
        );

        link_heat_transfer_entity(
            &mut front_cv,
            heat_trf_entity,
            heat_transfer_interaction,
        )?;

        // now that the front_cv is edited, we can extract it 

        let front_cv_edited: SingleCVNode = match front_cv {
            HeatTransferEntity::ControlVolume(cv_type) => {
                match cv_type {
                    CVType::SingleCV(single_cv) => {
                        single_cv
                    },
                    _ => unimplemented!(),
                }
            },
            _ => unimplemented!(),
        };

        self.front_single_cv = front_cv_edited;


        Ok(())
    }

    /// connects an adjacent heat transfer entity to front cv 
    ///
    pub fn link_heat_transfer_entity_to_back_cv(&mut self,
    heat_trf_entity: &mut HeatTransferEntity,
    heat_transfer_interaction: HeatTransferInteractionType) 
        -> Result<(), ThermalHydraulicsLibError>{

        // let's find the back cv, 
        // we shall clone it, link it, then replace the back cv 
        // with the edited control volume 

        let back_cv_clone = self.back_single_cv.clone();

        // package it as a heat transfer entity 

        let mut back_cv: HeatTransferEntity = 
        HeatTransferEntity::ControlVolume(
            CVType::SingleCV(back_cv_clone)
        );

        link_heat_transfer_entity(
            &mut back_cv,
            heat_trf_entity,
            heat_transfer_interaction,
        )?;

        // now that the back_cv is edited, we can extract it 

        let back_cv_edited: SingleCVNode = match back_cv {
            HeatTransferEntity::ControlVolume(cv_type) => {
                match cv_type {
                    CVType::SingleCV(single_cv) => {
                        single_cv
                    },
                    _ => unimplemented!(),
                }
            },
            _ => unimplemented!(),
        };

        self.back_single_cv = back_cv_edited;


        Ok(())
    }
}
