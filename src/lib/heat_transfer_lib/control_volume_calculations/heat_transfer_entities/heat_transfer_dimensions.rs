use uom::si::f64::*;
/// XThicknessThermalConduction is essentially a struct containing 
/// one length describing a thickness in cartesian coordinates 
/// for thermal conduction.
///
/// This type is meant for an input for various enums and functions 
/// in this crate.
///
/// It is meant to guide the user so that they know what the 
/// length inputs represents.
///
/// ```rust
///
/// use uom::si::length::meter;
/// use uom::si::f64::*;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// control_volume_calculations:: heat_transfer_entities:: 
/// XThicknessThermalConduction;
///
/// // let's say you have a thickness of 0.5 which you want to describe
///
/// let thickness_of_wall = Length::new::<meter>(0.5);
///
/// // we first need to convert it into an XThicknessThermalConduction
/// // type first
/// let wall_thickness_input = XThicknessThermalConduction::from(
/// thickness_of_wall);
///
/// // to convert it back into a length type, we use the into() method 
///
/// let thickness_wall_for_calculation: Length = 
/// wall_thickness_input.into();
/// 
/// // both these are the same
/// assert_eq!(thickness_of_wall, thickness_wall_for_calculation);
/// ```
/// 
/// 
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct XThicknessThermalConduction {
    thickness: Length,
}

impl From<Length> for XThicknessThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for XThicknessThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
}

/// This represents a thickness for radial conduction for 
/// cylindrical shell layers
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct RadialCylindricalThicknessThermalConduction {
    thickness: Length,
}

impl From<Length> for RadialCylindricalThicknessThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for RadialCylindricalThicknessThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
}

/// This represents an inner diameter  for radial conduction 
/// for spherical and cylindrical shell 
/// layers
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct InnerDiameterThermalConduction {
    thickness: Length,
}

impl From<Length> for InnerDiameterThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for InnerDiameterThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
}

/// This represents an outer diameter 
/// for radial conduction for spherical and Cylindrical shell 
/// layers
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct OuterDiameterThermalConduction {
    thickness: Length,
}

impl From<Length> for OuterDiameterThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for OuterDiameterThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
}

/// This represents an tube length 
/// ie. axial length for a cylindrical body
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CylinderLengthThermalConduction {
    thickness: Length,
}

impl From<Length> for CylinderLengthThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for CylinderLengthThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
}
