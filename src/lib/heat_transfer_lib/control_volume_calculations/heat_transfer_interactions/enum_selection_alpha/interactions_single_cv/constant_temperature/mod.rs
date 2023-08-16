use uom::si::f64::*;



use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;

mod advection;
use advection::*;

mod conductance;
use conductance::*;

#[inline]
pub fn calculate_single_cv_node_front_constant_temperature_back(
    boundary_condition_temperature: ThermodynamicTemperature,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<(), String> {
    // this code is pretty crappy but I'll match advection first

    match interaction {
        HeatTransferInteractionType::Advection(
        advection_dataset) => {

                // I'm mapping my own error to string, so off
                calculate_cv_front_bc_back_advection(
                    boundary_condition_temperature,
                    control_vol,
                    advection_dataset).map_err(
                        |error|{
                            error.to_string()
                        })?;
                return Ok(());
            },
        _ => (),
    }

    // if anything else, use conductance

    calculate_single_cv_node_constant_temperature_conductance(
        boundary_condition_temperature,
        control_vol,
        interaction).map_err(
            |error|{
                error.to_string()
            })?;

    return Ok(());
}

/// for connecting a bc to cv where 
///
/// (cv) ---------- (constant temperature bc)
#[inline]
pub fn calculate_constant_temperature_front_single_cv_back(
    control_vol: &mut SingleCVNode,
    boundary_condition_temperature: ThermodynamicTemperature,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    match interaction {
        HeatTransferInteractionType::Advection(
        advection_dataset) => {

                // I'm mapping my own error to string, so off
                calculate_bc_front_cv_back_advection(
                    control_vol,
                    boundary_condition_temperature,
                    advection_dataset).map_err(
                        |error|{
                            error.to_string()
                        })?;
                return Ok(());
            },
        _ => (),
    }
    // if anything else, use conductance

    calculate_single_cv_node_constant_temperature_conductance(
        boundary_condition_temperature,
        control_vol,
        interaction).map_err(
            |error|{
                error.to_string()
            })?;
    return Ok(());
}

