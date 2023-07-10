
use uom::si::f64::*;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::heat_transfer_dimensions::SurfaceArea;
use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::
heat_transfer_dimensions::XThicknessThermalConduction;

use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::
heat_transfer_dimensions::CrossSectionalArea;

use crate::heat_transfer_lib::thermophysical_properties::Material;

/// here we have a struct for dual Cartesian Thermal conduction 
/// in three dimensions
/// on
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataDualCartesianThermalConductanceThreeDimension{

    pub material_1: Material,
    pub material_2: Material,
    pub xs_area: CrossSectionalArea,
    pub thickness_1: XThicknessThermalConduction,
    pub thickness_2: XThicknessThermalConduction,

}

/// here we have a struct for simple convection resistance
/// in three dimensions
/// on
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataUserSpecifiedConvectionResistance{
    pub surf_area: SurfaceArea,
    pub heat_transfer_coeff: HeatTransfer,

}
