use crate::{heat_transfer_lib::control_volume_calculations::heat_transfer_entities::{HeatTransferEntity, CVType}, thermal_hydraulics_error::ThermalHydraulicsLibError};

use super::SingleCVNode;

impl TryFrom<HeatTransferEntity> for SingleCVNode {

    fn try_from(hte: HeatTransferEntity) -> Result<Self, Self::Error> {
        match hte {
            HeatTransferEntity::ControlVolume(cv) => {
                match cv {
                    CVType::SingleCV(single_cv) => {
                        Ok(single_cv)
                    },
                    _ => Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity),
                }
            },
            _ => Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity),
        }
    }

    type Error = ThermalHydraulicsLibError;
}

impl Into<HeatTransferEntity> for SingleCVNode {
    fn into(self) -> HeatTransferEntity {
        HeatTransferEntity::ControlVolume(
            CVType::SingleCV(
                self
            ))

    }
}
