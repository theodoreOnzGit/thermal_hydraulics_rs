use enums_alpha::data_enum_structs::DataAdvection;
use uom::num_traits::Zero;
use uom::si::f64::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
use crate::heat_transfer_lib::thermophysical_properties::density::try_get_rho;
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::try_get_h;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions
::enum_selection_alpha::interactions_single_cv::constant_temperature)
fn calculate_cv_front_bc_back_advection(
    boundary_condition_temperature: ThermodynamicTemperature,
    control_vol: &mut SingleCVNode,
    advection_data: DataAdvection
) -> Result<(), ThermalHydraulicsLibError>{

    let mass_flow_from_bc_to_cv = advection_data.mass_flowrate;

    let specific_enthalpy_cv: AvailableEnergy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    // for the constant temperature bc, 
    // I'll just assume the fluid flowing from the bc is the same as the 
    // fluid in the cv
    //
    // and the pressure is same as the cv

    let control_vol_material = control_vol.material_control_volume;
    let control_vol_pressure = control_vol.pressure_control_volume;

    let specific_enthalpy_bc: AvailableEnergy = try_get_h(
        control_vol_material,
        boundary_condition_temperature,
        control_vol_pressure
    )?;
    // calculate heat rate 

    let heat_flowrate_from_bc_to_cv: Power 
    = advection_heat_rate(mass_flow_from_bc_to_cv,
        specific_enthalpy_bc,
        specific_enthalpy_cv,)?;

    // push to cv
    control_vol.rate_enthalpy_change_vector.
        push(heat_flowrate_from_bc_to_cv);


    let density_cv = advection_data.fluid_density_heat_transfer_entity_2;

    // we also need to calculate bc density 
    // now this doesn't quite work well for compressible flow but for 
    // now I'll let it be

    let density_bc: MassDensity = try_get_rho(
        control_vol_material,
        boundary_condition_temperature,
        control_vol_pressure
    )?;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_bc_to_cv > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from bc 
        // to cv 
        // then the density we use is bc 

        volumetric_flowrate = mass_flow_from_bc_to_cv/density_bc;

    } else {
        // if mass flowrate is positive, flow is moving from cv
        // to bc
        // then the density we use is cv

        volumetric_flowrate = mass_flow_from_bc_to_cv/density_cv;
    }


    // for courant number
    control_vol.volumetric_flowrate_vector.push(
        volumetric_flowrate);


    Ok(())
}

#[inline]
pub (in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions
::enum_selection_alpha::interactions_single_cv::constant_temperature)
fn calculate_bc_front_cv_back_advection(
    control_vol: &mut SingleCVNode,
    boundary_condition_temperature: ThermodynamicTemperature,
    advection_data: DataAdvection
) -> Result<(), ThermalHydraulicsLibError>{

    let mass_flow_from_cv_to_bc = advection_data.mass_flowrate;

    let specific_enthalpy_cv: AvailableEnergy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    // for the constant temperature bc, 
    // I'll just assume the fluid flowing from the bc is the same as the 
    // fluid in the cv
    //
    // and the pressure is same as the cv

    let control_vol_material = control_vol.material_control_volume;
    let control_vol_pressure = control_vol.pressure_control_volume;

    let specific_enthalpy_bc: AvailableEnergy = try_get_h(
        control_vol_material,
        boundary_condition_temperature,
        control_vol_pressure
    )?;
    // calculate heat rate 

    let heat_flowrate_from_bc_to_cv: Power 
    = advection_heat_rate(mass_flow_from_cv_to_bc,
        specific_enthalpy_cv,
        specific_enthalpy_bc,)?;

    // push to cv
    control_vol.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_bc_to_cv);


    let density_cv = advection_data.fluid_density_heat_transfer_entity_2;

    // we also need to calculate bc density 
    // now this doesn't quite work well for compressible flow but for 
    // now I'll let it be

    let density_bc: MassDensity = try_get_rho(
        control_vol_material,
        boundary_condition_temperature,
        control_vol_pressure
    )?;

    let volumetric_flowrate: VolumeRate;

    if mass_flow_from_cv_to_bc > MassRate::zero() {
        // if mass flowrate is positive, flow is moving from bc 
        // to cv 
        // then the density we use is bc 

        volumetric_flowrate = mass_flow_from_cv_to_bc/density_bc;

    } else {
        // if mass flowrate is positive, flow is moving from cv
        // to bc
        // then the density we use is cv

        volumetric_flowrate = mass_flow_from_cv_to_bc/density_cv;
    }


    // for courant number
    control_vol.volumetric_flowrate_vector.push(
        -volumetric_flowrate);


    Ok(())
}
