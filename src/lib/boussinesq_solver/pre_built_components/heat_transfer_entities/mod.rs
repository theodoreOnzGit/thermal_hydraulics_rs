use self::cv_types::CVType;
use self::bc_types::BCType;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
/// Contains entities which transfer heat and interact with each 
/// other
///
/// for example, control volumes and boundary conditions
#[derive(Debug,Clone,PartialEq)]
pub enum HeatTransferEntity {
    /// Contains a list of ControlVolumeTypes
    ControlVolume(CVType),
    /// Contains a list of Boundary conditions
    BoundaryConditions(BCType)
}

impl HeatTransferEntity {

    /// allows the user to override the heat transfer entity 
    pub fn set(&mut self, 
        user_input_hte: HeatTransferEntity) -> Result<(), ThermalHydraulicsLibError>{
        *self = user_input_hte;

        Ok(())
    }
}

/// all the types of Control volumes are represented in an enum 
/// to abstract away the complications of connecting different types 
/// of control volumes. 
pub mod cv_types;

/// all the types of boundary conditions are represented in an enum 
/// to abstract away the complications of connecting different types 
/// of boundary conditions
pub mod bc_types;
