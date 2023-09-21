use super::data_enum_structs::DataAdvection;
use super::HeatTransferInteractionType;
use uom::si::f64::*;

impl HeatTransferInteractionType {

    /// constructs a new advection interaction so it's less 
    /// cumbersome for the user
    pub fn new_advection_interaction(
    mass_flowrate: MassRate,
    fluid_density_heat_transfer_entity_1: MassDensity,
    fluid_density_heat_transfer_entity_2: MassDensity) -> Self {

        let advection_interaction: HeatTransferInteractionType = 
        HeatTransferInteractionType::Advection(
            DataAdvection {
                mass_flowrate,
                fluid_density_heat_transfer_entity_1,
                fluid_density_heat_transfer_entity_2,
            }
        );

        return advection_interaction;

    }
}
