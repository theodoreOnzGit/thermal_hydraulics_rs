use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::CVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;

use super::FluidArray;

impl Into<ArrayCVType> for FluidArray {
    fn into(self) -> ArrayCVType {
        ArrayCVType::GenericPipe(self)
    }
}

impl Into<CVType> for FluidArray {
    fn into(self) -> CVType {
        CVType::ArrayCV(
            ArrayCVType::GenericPipe(
                self))
    }
}

impl Into<HeatTransferEntity> for FluidArray {
    fn into(self) -> HeatTransferEntity {
        HeatTransferEntity::ControlVolume(
            CVType::ArrayCV(
                ArrayCVType::GenericPipe(
                    self)))
    }
}
