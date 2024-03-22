use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::bc_types::BCType;
use super::cv_types::CVType;
use super::HeatTransferEntity;


impl Into<HeatTransferEntity> for BCType {
    fn into(self) -> HeatTransferEntity {
        HeatTransferEntity::BoundaryConditions(self)
    }
}


impl TryFrom<HeatTransferEntity> for BCType {
    type Error = ThermalHydraulicsLibError;

    fn try_from(hte: HeatTransferEntity) -> Result<Self, Self::Error> {
        let bc_type: BCType = match hte {
            HeatTransferEntity::ControlVolume(_) => {
                return Err(
                    ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity)
            },
            HeatTransferEntity::BoundaryConditions(bc_type) => {
                bc_type
            },
        };

        Ok(bc_type)
    }
}

impl Into<HeatTransferEntity> for CVType {
    fn into(self) -> HeatTransferEntity {
        HeatTransferEntity::ControlVolume(self)
    }
}

impl TryFrom<HeatTransferEntity> for CVType {
    type Error = ThermalHydraulicsLibError;

    fn try_from(value: HeatTransferEntity) -> Result<Self, Self::Error> {
        match value {
            HeatTransferEntity::ControlVolume(cv) => {
                return Ok(cv);
            },
            HeatTransferEntity::BoundaryConditions(_) => {
                return Err(ThermalHydraulicsLibError::TypeConversionErrorHeatTransferEntity);
            },
        }
    }
}
