use uom::si::f64::*;

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct InnerDiameterThermalConduction {
    od: Length,
}

impl From<Length> for InnerDiameterThermalConduction {
    fn from(od: Length) -> Self{
        Self { od }
    }
}

impl Into<Length> for InnerDiameterThermalConduction {
    fn into(self) -> Length {
        self.od
    }
}
/// This represents an outer diameter 
/// for radial conduction for spherical and Cylindrical shell 
/// layers
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct OuterDiameterThermalConduction {
    od: Length,
}

impl From<Length> for OuterDiameterThermalConduction {
    fn from(od: Length) -> Self{
        Self { od }
    }
}

impl Into<Length> for OuterDiameterThermalConduction {
    fn into(self) -> Length {
        self.od
    }
}

/// This represents an tube length 
/// ie. axial length for a cylindrical body
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct CylinderLengthThermalConduction {
    cylinder_length: Length,
}

impl From<Length> for CylinderLengthThermalConduction {
    fn from(cylinder_length: Length) -> Self{
        Self { cylinder_length }
    }
}

impl Into<Length> for CylinderLengthThermalConduction {
    fn into(self) -> Length {
        self.cylinder_length
    }
}
/// for 1D calcs, we use a unit area
pub const UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS: f64 = 1.0;
