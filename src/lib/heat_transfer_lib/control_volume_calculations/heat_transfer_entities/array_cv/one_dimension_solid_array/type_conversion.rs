use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::CVType;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::SolidColumn;

impl Into<ArrayCVType> for SolidColumn {
    fn into(self) -> ArrayCVType {
        ArrayCVType::GenericColumn(self)
    }
}

impl Into<CVType> for SolidColumn {
    fn into(self) -> CVType {
        CVType::ArrayCV(
            ArrayCVType::GenericColumn(
                self))
    }
}

impl Into<HeatTransferEntity> for SolidColumn {
    fn into(self) -> HeatTransferEntity {
        HeatTransferEntity::ControlVolume(
            CVType::ArrayCV(
                ArrayCVType::GenericColumn(
                    self)))
    }
}

impl TryFrom<ArrayCVType> for SolidColumn {
    type Error = ThermalHydraulicsLibError;

    fn try_from(array_cv: ArrayCVType) -> Result<Self, Self::Error> {
        match array_cv {
            ArrayCVType::GenericColumn(solid_cv) => {
                Ok(solid_cv)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }
}

impl TryFrom<CVType> for SolidColumn {
    type Error = ThermalHydraulicsLibError;

    fn try_from(cv_type: CVType) -> Result<Self, Self::Error> {
        match cv_type {
            CVType::ArrayCV(array_cv_type) => {
                let fluid_array: SolidColumn 
                = array_cv_type.try_into()?;
                Ok(fluid_array)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }
}

impl TryFrom<HeatTransferEntity> for SolidColumn {
    type Error = ThermalHydraulicsLibError;

    fn try_from(hte: HeatTransferEntity) -> Result<Self, Self::Error> {
        match hte {
            HeatTransferEntity::ControlVolume(cv) => {
                let cv: SolidColumn 
                =  cv.try_into()?;
                Ok(cv)
            },
            _ => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
        }
    }
}
