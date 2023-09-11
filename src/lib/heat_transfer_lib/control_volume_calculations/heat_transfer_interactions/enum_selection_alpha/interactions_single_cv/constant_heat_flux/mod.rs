use std::f64::consts::PI;
use uom::si::f64::*;



use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;

mod advection;
use advection::*;

/// calculates the interaction between a heat flux BC and 
/// a control volume 
///
/// (single cv) ------------------ (heat flux bc)
///
/// the heat addition is at the front, the cv is at the back
pub fn calculate_constant_heat_flux_front_single_cv_back(
    control_vol: &mut SingleCVNode,
    heat_flux_into_control_vol: HeatFluxDensity,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    // first, obtain a heat transfer area from the constant heat flux 
    // BC
    let heat_transfer_area: Area = match interaction{
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            DualCartesianThermalConductance(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(_, _, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or \n
                Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        // these interaction types are acceptable
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => area,

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(
            cylinder_length, od) => {
                let od: Length = od.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * od * cylinder_length;
                area
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(
            cylinder_length, id) => {
                let id: Length = id.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * id * cylinder_length;
                area

            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        HeatTransferInteractionType::
            Advection(advection_data) => {

                calculate_bc_front_cv_back_advection(
                    control_vol,
                    advection_data).map_err(
                        |err| err.to_string()
                    )?;

                return Ok(());
            },
    };

    let heat_flowrate_into_control_vol: Power = 
    heat_flux_into_control_vol * heat_transfer_area;

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_into_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change

    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    //
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;
            ()
        },
        Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }
    
    return Ok(());


}

/// calculates the interaction between a heat flux BC and 
/// a control volume 
///
/// (heat flux bc) ------------------ (single cv)
///
/// the cv is at the front 
/// heat addition is at the back
pub fn calculate_single_cv_front_heat_flux_back(
    heat_flux_into_control_vol: HeatFluxDensity,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    // first, obtain a heat transfer area from the constant heat flux 
    // BC
    let heat_transfer_area: Area = match interaction{
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            DualCartesianThermalConductance(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(_, _, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or \n
                Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(_, _) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(_, _) => 
            return Err("please specify interaction \n 
                type as UserSpecifiedHeatFluxCustomArea \n 
                or Similar".to_string()),

        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        // these interaction types are acceptable
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => area,

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(
            cylinder_length, od) => {
                let od: Length = od.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * od * cylinder_length;
                area
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(
            cylinder_length, id) => {
                let id: Length = id.into();
                let cylinder_length: Length  = cylinder_length.into();

                let area = PI * id * cylinder_length;
                area

            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => 
            return Err("please specify interaction type as \n 
                UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        HeatTransferInteractionType::
            Advection(advection_data) => 
            {
                calculate_cv_front_bc_back_advection(
                    control_vol,
                    advection_data).map_err(
                        |err| err.to_string()
                    )?;

                return Ok(());
            },
    };

    let heat_flowrate_into_control_vol: Power = 
    heat_flux_into_control_vol * heat_transfer_area;

    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_into_control_vol);

    // auto time stepping doesn't work for constant heat flux 
    // or specified power as well. 
    // it is best to see at the end of all power calculations what 
    // is the temperature change

    // For liquid CV, still need to calculate time scale based 
    // on convection flow
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
    //
    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    let cv_material = control_vol.material_control_volume;
    match cv_material {
        Solid(_) => {
            // in this case, we just have one cv and one bc 
            // so we only consider thermal inertia of this cv 
            calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                control_vol,
                interaction)?;
            ()
        },
        Liquid(_) => {
            // liquid time scales should be calculated using courant 
            // number at the end of each timestep after volumetric flows 
            // in and out of the cv are calculated
            ()
        },
    }
    
    return Ok(());


}
