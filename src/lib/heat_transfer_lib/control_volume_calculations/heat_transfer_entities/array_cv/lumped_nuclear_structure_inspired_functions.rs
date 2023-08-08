use ndarray::*;
use ndarray_linalg::*;
use uom::si::f64::*;
use uom::si::mass::kilogram;
use uom::si::power::watt;
use uom::si::thermal_conductance::watt_per_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::time::second;
use uom::ConstZero;
use uom::si::volume::cubic_meter;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::SingleCVNode;
use crate::heat_transfer_lib::thermophysical_properties::Material;
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::specific_enthalpy;

use super::solve_conductance_matrix_power_vector;
/// This is mostly a direct translation of GeN-Foam code, 
///
/// Now, this code is meant to calculate the internal temperature 
/// profile given an innermost temperature node, as well as
/// an outer cooling boundary condition
///
/// Here is the diagram provided from lumpedNuclearStructure.H:
/// Tmax            T[0]           T[1]          T[n-1]         Tsurface      
///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
///
///
/// for this code, i intend to adapt it as follows
///
/// Tmax            T[0]           T[1]          T[n-1]         T_ambient
/// [ignored]       T_singleCV                              
///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
///        [ignored]
///
/// There should be SingleCVNode on the inside, which gives the array 
/// cv a way to interact with other single CV nodes or other 
/// HeatTransferEntity objects
///
/// you must feed in a arrays of nodesNumber length long 
/// and nodesNumber + 1 for the thermal conductance arrays
#[inline]
#[allow(non_snake_case)]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn advance_timestep_for_externally_cooled_array_cv_no_insulation(
    inner_single_cv: &mut SingleCVNode,
    nodesNumber: usize,
    dt: Time,
    total_volume: Volume,
    q: Power,
    Tcool: ThermodynamicTemperature,
    Hcool: ThermalConductance,
    TOld: &mut Array1<ThermodynamicTemperature>,
    Hs: &mut Array1<ThermalConductance>, // nodesNumber +1 elements
    volFraction: &mut Array1<f64>,
    rhoCp: &mut Array1<VolumetricHeatCapacity>,
    qFraction: &mut Array1<f64>,
) 
-> Result<Array1<ThermodynamicTemperature>,error::LinalgError>{

    // First things first, we need to set up 
    // how the CV interacts with the internal array
    // here is heat added to CV
    

    let rate_enthalpy_change_vector: Vec<Power> = 
    inner_single_cv.rate_enthalpy_change_vector.clone();
    
    // I'm going to compute the power source for the inner most node

    let mut total_enthalpy_rate_change_innermost_node = 
    Power::new::<watt>(0.0);

    for enthalpy_chg_rate in 
        rate_enthalpy_change_vector.iter() {

            total_enthalpy_rate_change_innermost_node += *enthalpy_chg_rate;
        }

    // this innermost node will be an extra term added to the 
    // heat source vector S
    //
    // The new temperature will need to be used to calculate 
    // new specific enthalpy for the system

    let mut T: Array1<ThermodynamicTemperature> = TOld.map (
        |temp_kelvin_ptr: &ThermodynamicTemperature| {

            return *temp_kelvin_ptr;
        }

    );

    

    

    // now we start calculation

    if nodesNumber > 1 {
        let mut M: Array2<ThermalConductance> = 
        Array::zeros((nodesNumber, nodesNumber));
        let mut S: Array1<Power> = Array::zeros(nodesNumber);

        //
        // now, Hs is the conductance vector, but in GeN-Foam this is done 
        // on a per unit volume basis, I'm not quite going to do that 
        //
        // second note, the time scheme for this is implicit, not explicit
        // this makes this scheme quite stable, and larger fourier numbers 
        // are allowed for timestepping, good news!
        //
        // However, one must note that heat capacitance and conductance 
        // are not changing during the timestep, they are based on the 
        // old temperatures, so this doesn't quite work as well 
        // max mesh fourier number should still be 0.25
        //
        // note that Hs[0] is never used, it may as well be 0
        {
            M[[0,1]] = -Hs[1];
            M[[0,0]] = volFraction[0] * rhoCp[0] * total_volume 
                / dt + Hs[1];
            S[0] = q * qFraction[0] + 
                total_enthalpy_rate_change_innermost_node +
                TOld[0] * total_volume *
                volFraction[0] * rhoCp[0] /dt;
        }
        // Bulk medium
        if nodesNumber > 2 {
            for i in 1..nodesNumber-1 {

                M[[i,i+1]] = -Hs[i+1];
                M[[i,i-1]] = -Hs[i];
                M[[i,i]] = volFraction[i] * rhoCp[i] * total_volume / 
                    dt + Hs[i+1] + Hs[i];
                S[i] = q * qFraction[i] + TOld[i] * volFraction[i] * 
                    total_volume * rhoCp[i] / dt;
            }
        }

        {
            let i = nodesNumber-1;
            // this represents the total conductance to the outer 
            // node
            let HtoCool: ThermalConductance = 
            Hs[i+1]*Hcool/(Hs[i+1]+Hcool);
            M[[i,i-1]] = volFraction[i] * rhoCp[i] * total_volume / 
                dt + Hs[i] + HtoCool;
            M[[i,i]] = volFraction[i] * total_volume * rhoCp[i] / dt 
                + Hs[i] + HtoCool;
            S[i] = q * qFraction[i] 
                + TOld[i] * volFraction[i] * rhoCp[i] * total_volume / dt 
                + HtoCool * Tcool;       
        }

        T = solve_conductance_matrix_power_vector(M,S)?;


    } else {
        let HtoCool: ThermalConductance = Hs[1]*Hcool/(Hs[1]+Hcool);
        let M = volFraction[0] * total_volume * rhoCp[0] / dt + HtoCool;
        let S = q * qFraction[0] 
        + TOld[0] * volFraction[0] * rhoCp[0] * total_volume / dt 
        + HtoCool * Tcool;

        let temperature: TemperatureInterval = S/M;
        T[0] = ThermodynamicTemperature::new::<kelvin>( 
            temperature.value);
    }

    // Todo: probably need to synchronise error types in future
    let inner_node_enthalpy_next_timestep: AvailableEnergy = 
    specific_enthalpy(
        inner_single_cv.material_control_volume,
        T[0],
        inner_single_cv.pressure_control_volume).unwrap();

    inner_single_cv.current_timestep_control_volume_specific_enthalpy 
    = inner_node_enthalpy_next_timestep;

    // I also need to update the TOld vector 
    // This will ensure that the current temperature of the single 
    // cv node is equal to that of the matrix

    *TOld = T.mapv(
        |temperature_value| {
            return temperature_value;
        }
    );

    // set liquid cv mass 
    // probably also need to update error types in future
    inner_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
    inner_single_cv.rate_enthalpy_change_vector.clear();
    inner_single_cv.max_timestep_vector.clear();


    return Ok(T);
    
}


/// This is mostly a direct translation of GeN-Foam code, 
///
/// Now, this code is meant to calculate the internal temperature 
/// profile given an innermost temperature node, as well as
/// an outer cooling boundary condition
///
/// Here is the diagram provided from lumpedNuclearStructure.H:
/// Tmax            T[0]           T[1]          T[n-1]         Tsurface      
///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
///
///
/// for this code, i intend to adapt it as follows
///
/// Tmax            T[0]           T[1]          T[n-1]         T_ambient
/// [ignored]       T_innersingleCV             T_outersingleCV [ignored]
///   | --- H[0] --- | --- H[1] --- | --- H[2] --- |  --- H[n] --- |  
///        [ignored]                                    [ignore]
///
/// There should be SingleCVNode on the inside, which gives the array 
/// cv a way to interact with other single CV nodes or other 
/// HeatTransferEntity objects
///
/// you must feed in a arrays of nodesNumber length long 
/// and nodesNumber + 1 for the thermal conductance arrays
///
/// for this to work, we must have 2 nodes at least
///
/// What this means is that Hs[0] and Hs[n] both are zero
/// we don't even calculate T_ambient
#[inline]
#[allow(non_snake_case)]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn advance_timestep_for_conductance_array_cv(
    inner_single_cv: &mut SingleCVNode,
    outer_single_cv: &mut SingleCVNode,
    nodesNumber: usize,
    dt: Time,
    total_volume: Volume,
    q: Power,
    TOld: &mut Array1<ThermodynamicTemperature>,
    Hs: &mut Array1<ThermalConductance>, // nodesNumber +1 elements
    volFraction: &mut Array1<f64>,
    rhoCp: &mut Array1<VolumetricHeatCapacity>,
    qFraction: &mut Array1<f64>,
) 
-> Result<Array1<ThermodynamicTemperature>,error::LinalgError>{

    // First things first, we need to set up 
    // how the CV interacts with the internal array
    // here is heat added to CV
    

    let rate_enthalpy_change_vector_inner_node: Vec<Power> = 
    inner_single_cv.rate_enthalpy_change_vector.clone();
    
    // I'm going to compute the power source for the inner most node

    let mut total_enthalpy_rate_change_innermost_node = 
    Power::new::<watt>(0.0);

    for enthalpy_chg_rate in 
        rate_enthalpy_change_vector_inner_node.iter() {

            total_enthalpy_rate_change_innermost_node += *enthalpy_chg_rate;
        }

    let rate_enthalpy_change_vector_inner_node: Vec<Power> = 
    inner_single_cv.rate_enthalpy_change_vector.clone();
    
    // do the same for the outermost node
    let rate_enthalpy_change_vector_outer_node: Vec<Power> = 
    outer_single_cv.rate_enthalpy_change_vector.clone();
    
    // I'm going to compute the power source for the outer most node

    let mut total_enthalpy_rate_change_outermost_node = 
    Power::new::<watt>(0.0);

    for enthalpy_chg_rate in 
        rate_enthalpy_change_vector_outer_node.iter() {

            total_enthalpy_rate_change_outermost_node += *enthalpy_chg_rate;
        }

    let rate_enthalpy_change_vector_outer_node: Vec<Power> = 
    outer_single_cv.rate_enthalpy_change_vector.clone();

    // stopped here, need a break


    let mut T: Array1<ThermodynamicTemperature> = TOld.map (
        |temp_kelvin_ptr: &ThermodynamicTemperature| {

            return *temp_kelvin_ptr;
        }

    );

    

    

    // now we start calculation

    if nodesNumber > 1 {
        let mut M: Array2<ThermalConductance> = 
        Array::zeros((nodesNumber, nodesNumber));
        let mut S: Array1<Power> = Array::zeros(nodesNumber);

        //
        // now, Hs is the conductance vector, but in GeN-Foam this is done 
        // on a per unit volume basis, I'm not quite going to do that 
        //
        // second note, the time scheme for this is implicit, not explicit
        // this makes this scheme quite stable, and larger fourier numbers 
        // are allowed for timestepping, good news!
        //
        // However, one must note that heat capacitance and conductance 
        // are not changing during the timestep, they are based on the 
        // old temperatures, so this doesn't quite work as well 
        // max mesh fourier number should still be 0.25
        //
        // note that Hs[0] is never used, it may as well be 0
        {
            M[[0,1]] = -Hs[1];
            M[[0,0]] = volFraction[0] * rhoCp[0] * total_volume 
                / dt + Hs[1];
            S[0] = q * qFraction[0] + 
                total_enthalpy_rate_change_innermost_node +
                TOld[0] * total_volume *
                volFraction[0] * rhoCp[0] /dt;
        }
        // Bulk medium
        if nodesNumber > 2 {
            for i in 1..nodesNumber-1 {

                M[[i,i+1]] = -Hs[i+1];
                M[[i,i-1]] = -Hs[i];
                M[[i,i]] = volFraction[i] * rhoCp[i] * total_volume / 
                    dt + Hs[i+1] + Hs[i];
                S[i] = q * qFraction[i] + TOld[i] * volFraction[i] * 
                    total_volume * rhoCp[i] / dt;
            }
        }

        {
            let i = nodesNumber-1;
            // this represents the total conductance to the outer 
            // node
            //  at the ambient temperature
            //
            //  so i'm not going to edit the equations
            //  but I set HtoCool to zero 
            //  and Tcool to room temp 
            //  doesn't matter what Tcool is because 
            //  the terms become zero anyhow
            let HtoCool: ThermalConductance = 
            ThermalConductance::new::<watt_per_kelvin>(0.0);
            let Tcool = ThermodynamicTemperature::new::<kelvin>(293.0);

            M[[i,i-1]] = volFraction[i] * rhoCp[i] * total_volume / 
                dt + Hs[i] + HtoCool;
            M[[i,i]] = volFraction[i] * total_volume * rhoCp[i] / dt 
                + Hs[i] + HtoCool;
            S[i] = q * qFraction[i] 
                + TOld[i] * volFraction[i] * rhoCp[i] * total_volume / dt 
                + HtoCool * Tcool;       
        }

        T = solve_conductance_matrix_power_vector(M,S)?;


    } else {
        panic!("nodesNumber must be at least 2");

        // legacy code 
        //
        //let HtoCool: ThermalConductance = Hs[1]*Hcool/(Hs[1]+Hcool);
        //let M = volFraction[0] * total_volume * rhoCp[0] / dt + HtoCool;
        //let S = q * qFraction[0] 
        //+ TOld[0] * volFraction[0] * rhoCp[0] * total_volume / dt 
        //+ HtoCool * Tcool;

        //let temperature: TemperatureInterval = S/M;
        //T[0] = ThermodynamicTemperature::new::<kelvin>( 
        //    temperature.value);
    }

    // Todo: probably need to synchronise error types in future
    let inner_node_enthalpy_next_timestep: AvailableEnergy = 
    specific_enthalpy(
        inner_single_cv.material_control_volume,
        T[0],
        inner_single_cv.pressure_control_volume).unwrap();

    inner_single_cv.current_timestep_control_volume_specific_enthalpy 
    = inner_node_enthalpy_next_timestep;

    // I also need to update the TOld vector 
    // This will ensure that the current temperature of the single 
    // cv node is equal to that of the matrix

    *TOld = T.mapv(
        |temperature_value| {
            return temperature_value;
        }
    );

    // set liquid cv mass 
    // probably also need to update error types in future
    inner_single_cv.set_liquid_cv_mass_from_temperature().unwrap();
    inner_single_cv.rate_enthalpy_change_vector.clear();
    inner_single_cv.max_timestep_vector.clear();


    return Ok(T);
    
}
