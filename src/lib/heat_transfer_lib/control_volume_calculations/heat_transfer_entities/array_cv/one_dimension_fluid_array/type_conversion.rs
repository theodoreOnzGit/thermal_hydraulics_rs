use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::CVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

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

impl TryFrom<ArrayCVType> for FluidArray {
    type Error = ThermalHydraulicsLibError;

    fn try_from(array_cv: ArrayCVType) -> Result<Self, Self::Error> {
        match array_cv {
            ArrayCVType::GenericPipe(fluid_cv) => {
                Ok(fluid_cv)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }
}

impl TryFrom<CVType> for FluidArray {
    type Error = ThermalHydraulicsLibError;

    fn try_from(cv_type: CVType) -> Result<Self, Self::Error> {
        match cv_type {
            CVType::ArrayCV(array_cv_type) => {
                let fluid_array: FluidArray 
                = array_cv_type.try_into()?;
                Ok(fluid_array)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }
}

impl TryFrom<HeatTransferEntity> for FluidArray {
    type Error = ThermalHydraulicsLibError;

    fn try_from(hte: HeatTransferEntity) -> Result<Self, Self::Error> {
        match hte {
            HeatTransferEntity::ControlVolume(cv) => {
                let cv: FluidArray 
                =  cv.try_into()?;
                Ok(cv)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }

}

