use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib;



use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
/// calculates a conductance interaction between the constant 
/// temperature bc and cv
///
/// for conductance, orientation of bc and cv does not usually matter
#[inline]
pub(in crate) 
fn calculate_single_cv_node_constant_temperature_conductance(
    boundary_condition_temperature: ThermodynamicTemperature,
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError> {
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



    
    // we'll need thermal conductance otherwise 
    let cv_bc_conductance: ThermalConductance = 
    get_thermal_conductance_based_on_interaction(
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
