use uom::si::f64::*;
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

impl CVType {
    #[inline]
    /// gets the material 
    pub fn get_material(&mut self) -> Result<Material,ThermalHydraulicsLibError>{


        match self {
            CVType::SingleCV(single_cv_node) => {
                return Ok(single_cv_node.material_control_volume);
            },
            CVType::ArrayCV(array_cv) => {
                array_cv.get_array_cv_material()
            },
        }
    }

    #[inline]
    fn get_temperature_vector(&mut self) -> 
    Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        match self {
            CVType::SingleCV(single_cv) => {
                let temperature = single_cv.get_temperature_from_enthalpy_and_set()?;

                let mut temp_vec: Vec<ThermodynamicTemperature> = vec![];

                temp_vec.push(temperature);

                return Ok(temp_vec);

            },
            CVType::ArrayCV(array_cv) => {
                return array_cv.get_temperature_vector();
            },
        }
    }
}



/// contains codes for boundary conditions
pub mod boundary_conditions;
pub use boundary_conditions::BCType;



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
pub use array_cv::*;


/// preprocessing contains auto timestepping and other functions
pub mod preprocessing;


/// postprocessing contains functions to obtain temperature profiles 
/// of the HeatTransferEntity
pub mod postprocessing;


/// calculation modules contain methods to advance timestep 
pub mod calculation;


use crate::heat_transfer_lib::thermophysical_properties::Material;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;



