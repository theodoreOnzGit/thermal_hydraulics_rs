
use std::f64::consts::PI;

use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;

use super::heat_transfer_entities::HeatTransferEntity;
use uom::si::f64::*;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::pressure::atmosphere;
/// represents heater version 2 without insulation 
/// This is because during 2018-ish, the heater insulation 
/// got burnt off and a lot of frequency response tests were done 
/// with insulation removed
///
/// Heater version 2 bare has no insulation
/// but it has a twisted tape interior
///
///
/// note that it only contains the heated section, not the top nor 
/// bottom heads
///
/// Note: need to check for memory leaks
#[derive(Debug,Clone,PartialEq)]
pub struct StructuralSupport {

    inner_nodes: usize,

    pub support_array: HeatTransferEntity,

    pub ambient_temperature: ThermodynamicTemperature,

    pub heat_transfer_to_air: HeatTransfer,

    pub total_lateral_surface_area: Area,


}

impl StructuralSupport {


    /// constructs a structural support made typically of steel
    /// shaped in a cylinder
    /// for simplicity
    ///
    /// Unheated Structure Thermal Inertia: ignored
    pub fn new_steel_support_cylinder(
        component_length: Length,
        diameter: Length,
        initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature) -> Self {

        let user_specified_inner_nodes: usize = 0;
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let total_lateral_surface_area = 
        diameter * PI * component_length;

        // theoretically it's 6 W/(m^2 K) but then we'll have to manually 
        // input wall structures for additional heat loss
        //
        let h_to_air: HeatTransfer = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(6.0);

        // correlation 




        // now the outer steel array
        let steel_shell_array = 
        SolidColumn::new_cylinder(
            component_length,
            diameter,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { 
            inner_nodes: user_specified_inner_nodes,
            support_array: steel_shell_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
            total_lateral_surface_area,
        };
    }

}




pub mod preprocessing;

pub mod calculation;

pub mod postprocessing;
