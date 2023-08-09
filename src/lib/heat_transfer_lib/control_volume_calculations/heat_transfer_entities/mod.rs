use std::f64::consts::PI;

use uom::si::{f64::*, pressure::atmosphere, power::watt, time::second, length::meter};

use crate::heat_transfer_lib::
thermophysical_properties::{Material, 
    specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy}, density::density, thermal_diffusivity::thermal_diffusivity, specific_heat_capacity::specific_heat_capacity};

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
/// placeholder, I'd like to have some associated functions to 
/// deal with the HeatTransferEntity type
///
/// probably one to get the courant number, 
/// and second, to use a timestep to calculate the new enthalpy 
/// and update enthalpy
/// 
/// last but not least, extract temperatures for sensing purposes
impl HeatTransferEntity {

    /// for control volumes, this method allows you to 
    /// calculate the enthalpy of the next timestep and 
    /// set the 
    /// current timestep enthalpy as the enthalpy calculated  
    /// for the next timestep
    ///
    /// you are required to explicitly provide a timestep for this 
    pub fn advance_timestep(entity: &mut HeatTransferEntity,
    timestep: Time) -> Result<(), String> {

        // first match CV or BC, 
        // Boundary conditions don't need to advance timestep
        // so we can leave them be (it should return an Ok(()) value 
        // rather than an Err() value)

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => return Ok(()),
        };

        // once I have the cv_type enum, match it again

        let cv_advance_result: Result<(), String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.advance_timestep(timestep)
                },
                CVType::ArrayCV(_) => return Err("not implemented".to_string()),
                //CVType::CartesianConduction1DArray(cv) => {
                //    cv.advance_timestep(timestep)
                //},
            };

        return cv_advance_result;
    }

    /// gets the temperature of the HeatTransferEntity 
    /// usually control volume at the current timestep
    pub fn temperature(entity: &mut HeatTransferEntity) -> 
    Result<ThermodynamicTemperature, String> {

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => 
                return Err("getting temperature not \n 
                    implemented for BoundaryConditions".to_string()),
        };

        // once I have the cv_type enum, match it again

        let cv_temperature_result: 
        Result<ThermodynamicTemperature, String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.get_temperature()
                },
                CVType::ArrayCV(_) => return Err("not implemented".to_string()),
                //CVType::CartesianConduction1DArray(cv) => {
                //    cv.get_bulk_temperature()
                //},
            };

        return cv_temperature_result;
    }

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
                CVType::ArrayCV(_) => return Err("not implemented".to_string()),
            };

        return cv_timestep_result;
    }

}

/// To determine heat transfer between two control volumes or 
/// generally, two heat transfer entities, one must determine 
/// which control volume is in front and which is at the back
/// 
/// This type would tell the solver that this control volume is 
/// in the front
/// 
#[derive(Debug,Clone,PartialEq)]
pub struct FrontHeatTransferEntity {
    entity: HeatTransferEntity,
}

impl From<HeatTransferEntity> for FrontHeatTransferEntity {
    fn from(entity: HeatTransferEntity) -> Self{
        Self { entity }
    }
}

impl Into<HeatTransferEntity> for FrontHeatTransferEntity {
    fn into(self) -> HeatTransferEntity {
        self.entity
    }
}

/// To determine heat transfer between two control volumes or 
/// generally, two heat transfer entities, one must determine 
/// which control volume is in front and which is at the back
/// 
/// This type would tell the solver that this control volume is 
/// in the back
/// 
#[derive(Debug,Clone,PartialEq)]
pub struct BackHeatTransferEntity {
    entity: HeatTransferEntity,
}

impl From<HeatTransferEntity> for BackHeatTransferEntity {
    fn from(entity: HeatTransferEntity) -> Self{
        Self { entity }
    }
}

impl Into<HeatTransferEntity> for BackHeatTransferEntity {
    fn into(self) -> HeatTransferEntity {
        self.entity
    }
}

/// Contains Types of Control Volumes (CVs)
#[derive(Debug,Clone,PartialEq)]
pub enum CVType {
    /// This CV is the most basic,  it can be represented by a single 
    /// point or node
    SingleCV(SingleCVNode),
    /// Array CVs are collections of SingleCVs, 
    /// or discretised singleCVs
    /// but do not require the 
    /// user to manually specify the connections between the SingleCVs
    ArrayCV(ArrayCVType),
}

/// Contains types of array control volumes
#[derive(Debug,Clone,PartialEq)]
pub enum ArrayCVType {
    /// This is one type of array CV which contains a 1D 
    /// conduction model with one material type
    /// standby for implementation
    Cartesian1D(CartesianConduction1DArray),
}

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



/// Contains different length types for use in defining interactions 
/// between heat transfer entities
///
/// This is to make it clear to the user exactly what 
/// kind of length we are specifying 
///
/// For example, to define an annular hollow cylinder, we need three 
/// lengths:
///
/// the zLength of the cylinder,
/// inner diameter 
/// outer diameter
///
/// each of these have type Length. The user would have to read 
/// the documentation as to what kind of length is being 
/// specified by the user 
///
/// Hence, I'm making some extra types to force the compiler to tell you 
/// (the user) what kind of length you need to specify
///
/// 
///
pub mod heat_transfer_dimensions;
pub use heat_transfer_dimensions::*;

/// contains functions and methods for the single control volumes
pub mod single_cv_node;
pub use single_cv_node::*;

/// contains functions and methods for array control volumes 
pub mod array_cv;
/// the array control volume (cv) module is basically to speed up creation 
/// process of control volumes which require somewhat of a one dimensional 
/// mesh or series of one dimensional meshes 
///
/// of course, you could link up several single_cv_node objects and 
/// use that as the array control volume. However, the creation process is 
/// quite tedious. The only advantage perhaps is flexibility and 
/// customisability 
///
/// Some of the code here will be a direct translation of GeN-Foam 
/// code licensed under GPL 3.0 into rs
///
/// Hence, this entire library is licensed under GPL 3.0
pub use array_cv::*;
