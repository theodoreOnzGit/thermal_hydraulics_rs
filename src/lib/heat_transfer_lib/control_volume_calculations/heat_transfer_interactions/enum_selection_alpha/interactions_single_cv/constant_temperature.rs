use enums_alpha::data_enum_structs::DataAdvection;
use uom::num_traits::Zero;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib;



use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
use crate::heat_transfer_lib::thermophysical_properties::density::density;
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::specific_enthalpy;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
pub fn calculate_single_cv_node_front_constant_temperature_back(
    control_vol: &mut SingleCVNode,
    boundary_condition_temperature: ThermodynamicTemperature,
    interaction: HeatTransferInteractionType) -> Result<(), String> {

    // first let's get the control volume temperatures out
    
    let cv_enthalpy = 
    control_vol.current_timestep_control_volume_specific_enthalpy;

    let cv_material = control_vol.material_control_volume;

    let cv_pressure = control_vol.pressure_control_volume;

    let cv_temperature = heat_transfer_lib::thermophysical_properties:: 
        specific_enthalpy::temperature_from_specific_enthalpy(
            cv_material, 
            cv_enthalpy, 
            cv_pressure)?;

    // for now we assume the boundary condition pressure is the same 
    // as the control volume pressure, because pressure is not 
    // specified or anything

    let bc_pressure = cv_pressure.clone();

    // this code is pretty crappy but I'll match advection first

    match interaction {
        HeatTransferInteractionType::Advection(
        advection_dataset) => {

                // I'm mapping my own error to string, so off
                calculate_cv_front_bc_back_advection(
                    boundary_condition_temperature,
                    control_vol,
                    advection_dataset
                )
                    .map_err(
                        |error|{
                            error.to_string()
                        })?;
                ()
            },
        _ => (),
    }


    
    // we'll need thermal conductance otherwise 
    let cv_bc_conductance: ThermalConductance = 
    get_thermal_conductance(
        cv_temperature, 
        boundary_condition_temperature, 
        cv_pressure, 
        bc_pressure, 
        interaction)?;

    // with conductance settled, we should be able to calculate a power 

    let bc_temp_minus_cv_temp_kelvin: f64 = 
        boundary_condition_temperature.get::<kelvin>() - 
        cv_temperature.get::<kelvin>();

    let bc_temp_mins_cv_temp: TemperatureInterval = 
    TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(
        bc_temp_minus_cv_temp_kelvin);

    // heat flow to the destination 
    // be is proportional to -(T_final - T_initial) 
    // or -(T_destination - T_source)
    //
    // so if the destination is the boundary condition,
    //
    let heat_flowrate_from_cv_to_bc: Power = 
    - cv_bc_conductance * bc_temp_mins_cv_temp;

    // now, push the power change or heat flowrate 
    // to the control volume 
    //

    control_vol.rate_enthalpy_change_vector.
        push(-heat_flowrate_from_cv_to_bc);

    // for constant temperature BC in interaction with CV,
    // we only need take into consideration the CV 
    // timescale, this is based on the fourier number or Nusselt 
    // number edited fourier number
    //
    // for solids mesh fourier number need only 
    // be done once, not every time 
    // an interaction is formed 
    //
    // probably the cell stability fourier number will be done in the 
    // constructor. however, with convection, the time scale must be 
    // recalculated at every time step. so it really depends whether 
    // it's solid or fluid control volume
    // 

    // Actually, for solid CV, I will also need to recalculate time scale 
    // based on the material thermal thermal_diffusivity
    // match statement is meant to tell that liquid CVs are not quite 
    // ready for use
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

#[inline]
pub (in crate)
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

    let specific_enthalpy_bc: AvailableEnergy = specific_enthalpy(
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

    let density_bc: MassDensity = density(
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
pub (in crate)
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

    let specific_enthalpy_bc: AvailableEnergy = specific_enthalpy(
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

    let density_bc: MassDensity = density(
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
