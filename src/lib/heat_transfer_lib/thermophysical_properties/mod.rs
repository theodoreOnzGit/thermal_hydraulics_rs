//! This module contains a library of liquid and solid 
//! thermophysical properties

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

/// Contains a selection of liquids with predefined material properties
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum LiquidMaterial {
    /// therminol VP1 
    TherminolVP1,
    /// DowthermA, using 
    DowthermA
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
pub mod error_types;

