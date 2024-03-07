//! This module contains a library of liquid and solid 
//! thermophysical properties

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

/// basically,
/// insert this enum into a thermophysical property function 
/// or something
/// then it will extract the 
/// thermophysical property for you in unit safe method
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum Material {
    /// Contains a list of selectable solids
    Solid(SolidMaterial),
    /// Contains a list of selectable liquids
    Liquid(LiquidMaterial)
}

impl Default for Material{
    fn default() -> Self {
        // defaults to copper 
        return Self::Solid(SolidMaterial::Copper);
    }
}

impl TryInto<SolidMaterial> for Material {
    type Error = ThermalHydraulicsLibError;

    fn try_into(self) -> Result<SolidMaterial, Self::Error> {
        match self {
            Material::Solid(material) => {
                Ok(material)
            },
            Material::Liquid(_) => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorMaterial)
            },
        }
    }
}

impl TryInto<LiquidMaterial> for Material {
    type Error = ThermalHydraulicsLibError;

    fn try_into(self) -> Result<LiquidMaterial, Self::Error> {
        match self {
            Material::Solid(_) => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorMaterial)
            },
            Material::Liquid(material) => {
                Ok(material)
            },
        }
    }
}

/// Contains a selection of solids with predefined material properties
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum SolidMaterial {
    /// stainless steel 304 L, 
    /// material properties from 
    /// Graves, R. S., Kollie, T. G., McElroy, D. L., 
    /// & Gilchrist, K. E. (1991). The thermal conductivity of 
    /// AISI 304L stainless steel. International journal of 
    /// thermophysics, 12, 409-415.
    SteelSS304L,
    /// Copper material
    Copper,
    /// Fiberglass material
    Fiberglass
}

impl Into<Material> for SolidMaterial {
    fn into(self) -> Material {
        Material::Solid(self)
    }
}

/// Contains a selection of liquids with predefined material properties
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum LiquidMaterial {
    /// therminol VP1 
    TherminolVP1,
    /// DowthermA, using 
    DowthermA
}

impl Into<Material> for LiquidMaterial {
    fn into(self) -> Material {
        Material::Liquid(self)
    }
}

/// Density calculation
pub mod density;

/// Thermal conductivity calculation 
pub mod thermal_conductivity;

/// SpecificHeatCapacity calculation 
pub mod specific_heat_capacity;

/// dynamic viscosity calculation 
pub mod dynamic_viscosity;

/// specific enthalpy calculation
pub mod specific_enthalpy;

/// (hopefully) useful information in case a function fails to behave 
/// properly 
///
/// probably going to deprecate or integrate with existing 
/// thermal_hydraulics_error
pub mod error_types;

/// thermal diffusivity 
pub mod thermal_diffusivity;

/// volumetric_heat_capacity 
pub mod volumetric_heat_capacity;

/// prandtl number 
pub mod prandtl;

/// surface roughness 
pub mod solid_material_surface_roughness;

