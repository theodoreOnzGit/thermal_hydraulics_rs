
use uom::si::f64::*;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::HeatTransferEntity;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::heat_transfer_dimensions::SurfaceArea;
use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::
heat_transfer_dimensions::XThicknessThermalConduction;

use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::
heat_transfer_dimensions::CrossSectionalArea;

use crate::heat_transfer_lib::thermophysical_properties::LiquidMaterial;
use crate::heat_transfer_lib::thermophysical_properties::Material;

use super::HeatTransferInteractionType;

/// here we have a struct for dual Cartesian Thermal conduction 
/// in three dimensions
/// on
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataDualCartesianThermalConductanceThreeDimension{
    /// material for first cv
    pub material_1: Material,
    /// material for second cv
    pub material_2: Material,
    /// cross sectional area at interface
    pub xs_area: CrossSectionalArea,
    /// thickness of first cv
    pub thickness_1: XThicknessThermalConduction,
    /// thickness of second cv
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

impl DataAdvection {

    /// constructs an advection interaction by specifying 
    /// a fluid material 
    /// temperature of the heat transfer entity 1
    /// and temperature of heat transfer entity 2
    ///
    ///
    /// (heat tranfer entity 1) ----mass flowrate --> (heat transfer entity 2)
    ///
    #[inline]
    pub fn new_from_temperature_and_liquid_material(
        user_input_mass_flowrate: MassRate,
        fluid_material: LiquidMaterial,
        temperature_1: ThermodynamicTemperature,
        temperature_2: ThermodynamicTemperature,
    ) -> Self {


        let density_1 = fluid_material.density(temperature_1).unwrap();
        let density_2 = fluid_material.density(temperature_2).unwrap();
        return Self {
            mass_flowrate: user_input_mass_flowrate,
            fluid_density_heat_transfer_entity_1: density_1,
            fluid_density_heat_transfer_entity_2: density_2,
        };
    }

    /// constructs an advection interaction by specifying 
    /// a fluid material 
    /// mutable reference of the heat transfer entity 1
    /// and mutable reference of heat transfer entity 2
    ///
    ///
    /// (heat tranfer entity 1) ----mass flowrate --> (heat transfer entity 2)
    ///
    #[inline] 
    pub fn new_from_heat_transfer_entity(
        user_input_mass_flowrate: MassRate,
        fluid_material: LiquidMaterial,
        hte_1: &mut HeatTransferEntity,
        hte_2: &mut HeatTransferEntity) -> Self {

        // if one of the heat transfer entities is a adiabatic 
        // bc, it will have no temperature, but 
        // i can use temperature of the other heat transfer 
        // entity
        //
        // perhaps in future I can handle the errors better 
        //
        // for now, no choice but to unwrap due to time 
        // constraints
        let temperature_1: ThermodynamicTemperature = 
        HeatTransferEntity::temperature(hte_1).unwrap();

        let temperature_2: ThermodynamicTemperature = 
        HeatTransferEntity::temperature(hte_2).unwrap();

        let density_1 = fluid_material.density(temperature_1).unwrap();
        let density_2 = fluid_material.density(temperature_2).unwrap();

        return Self {
            mass_flowrate: user_input_mass_flowrate,
            fluid_density_heat_transfer_entity_1: density_1,
            fluid_density_heat_transfer_entity_2: density_2,
        };

    }
}

impl Into<HeatTransferInteractionType> for DataAdvection {
    fn into(self) -> HeatTransferInteractionType {
        HeatTransferInteractionType::Advection(self)
    }
}
