
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

    /// surface area for heat convection
    pub surf_area: SurfaceArea,
    /// heat transfer coefficient in watts per meter per kelvin
    pub heat_transfer_coeff: HeatTransfer,

}

/// here we have a useful for necessary advection information 

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataAdvection{

    /// mass flowrate
    pub mass_flowrate: MassRate,
    /// fluid density of control volume on left
    ///
    /// which means when you link control volumes or boundary 
    /// link(cv1, cv2, interaction)
    ///
    /// the picture is like this 
    ///
    /// (cv1) ----> advection ---> (cv2)
    ///
    /// cv1 is the left control volume 
    /// cv2 is the right control volume
    ///
    /// now, the cv is not always a cv, it could be any heat
    /// transfer entity
    pub fluid_density_heat_transfer_entity_1: MassDensity,
    /// fluid density of control volume on left
    ///
    /// which means when you link control volumes or boundary 
    /// link(cv1, cv2, interaction)
    ///
    /// the picture is like this 
    ///
    /// (cv1) ----> advection ---> (cv2)
    ///
    /// cv1 is the left control volume 
    /// cv2 is the right control volume
    /// now, the cv is not always a cv, it could be any heat
    /// transfer entity
    pub fluid_density_heat_transfer_entity_2: MassDensity

}
