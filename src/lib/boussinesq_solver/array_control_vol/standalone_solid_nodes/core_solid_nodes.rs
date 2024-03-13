use ndarray::*;
use ndarray_linalg::{*, error::LinalgError};
use uom::num_traits::Zero;
use uom::si::f64::*;
use uom::si::power::watt;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::specific_enthalpy::try_get_h;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use crate::boussinesq_solver::single_control_vol::SingleCVNode;
use crate::boussinesq_solver::array_control_vol::standalone_fluid_nodes::solve_conductance_matrix_power_vector;

use super::shell_solid_node::advance_timestep_solid_cylindrical_shell_node_no_axial_conduction;


/// calculates solid cylindrical core without axial conduction 
/// algorithm is simple, call on the shell_solid_node calculation 
/// but make one side have zero heat conductance
pub fn advance_timestep_solid_cylindrical_core_node_no_axial_conduction(
    number_of_nodes: usize,
    dt: Time,
    total_volume: Volume,
    q: Power,
    last_timestep_temperature_adjacent_side: &mut Array1<ThermodynamicTemperature>,
    solid_adjacent_side_conductance_array: &mut Array1<ThermalConductance>,
    last_timestep_temperature_solid: &mut Array1<ThermodynamicTemperature>,
    volume_fraction_array: &mut Array1<f64>,
    rho_cp: &mut Array1<VolumetricHeatCapacity>,
    q_fraction: &mut Array1<f64>)
-> Result<Array1<ThermodynamicTemperature>,ThermalHydraulicsLibError>{

    // I'm just re-calling the 
    // advance_timestep_solid_cylindrical_shell_node_no_axial_conduction
    // but with one side being an adiabatic bc (no conduction and no 
    // heat transfer)
    

    let mut adiabatic_conductance_array: Array1<ThermalConductance> = 
    Array::zeros(number_of_nodes);


    let adiabatic_conductance_array_mut_ref: &mut Array1<ThermalConductance> 
    = &mut adiabatic_conductance_array;
    

    // the temperature can be anything 

    let mut arbitrary_temperature_conductance_array: 
    Array1<ThermodynamicTemperature> = Array::default(number_of_nodes);

    let arbitrary_temperature_array_mut_ref: &mut Array1<ThermodynamicTemperature> 
    = &mut arbitrary_temperature_conductance_array;


    advance_timestep_solid_cylindrical_shell_node_no_axial_conduction(
        number_of_nodes,
        dt,
        total_volume,
        q,
        last_timestep_temperature_adjacent_side,
        solid_adjacent_side_conductance_array,
        arbitrary_temperature_array_mut_ref,
        adiabatic_conductance_array_mut_ref,
        last_timestep_temperature_solid,
        volume_fraction_array,
        rho_cp,
        q_fraction
    )
}
