use crate::{heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity, thermal_hydraulics_error::ThermalHydraulicsLibError};

use super::BCType;

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
