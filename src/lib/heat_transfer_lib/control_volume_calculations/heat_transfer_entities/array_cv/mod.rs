// to start off array_cv, we take inspiration from GeN-Foam's 
// lumpedNuclearStructure code 
//
// lumpedNuclearStructure was meant to model a heat generating 
// pebble geometry because GeN-Foam used to only allow 
// for pin shaped geometry. To deal with this, we have several 
// control volumes or nodes with thermal conductances between 
// the nodes
//
//
//
// now, for matrix solution, i use intel-mkl-static in ndarray_linalg 
// library 
//
// this is because ndarray_linalg using intel-mkl-static is cross 
// platform, and it can be used for windows, macos and linux 
//
// secondly, the intel-mkl-static library compiles the lapack library 
// locally and links it statically, rather than at the system level 
// therefore, the user won't have to worry as much about system 
// dependencies which can cause some headache
//
// btw, in future implementations of thermal_hydraulics_rs, i might 
// want to copy and paste how ndarray-linalg constructs its error types 
// so that my results are similar



use std::f64::consts::PI;

use ndarray::*;
use ndarray_linalg::*;
use uom::ConstZero;
use uom::si::f64::*;
use uom::si::power::watt;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::thermal_conductance::watt_per_kelvin;

/// This module contains direct translations of GeN-Foam code 
/// into rust for the lumped nuclear structure
#[path = "./gen-foam-lumped-nuclear-structure.rs"]
pub mod gen_foam_lumped_nuclear_structure;
use gen_foam_lumped_nuclear_structure::solve_conductance_matrix_power_vector;
/// This module contains adapts the GeN-Foam code for the 
/// lumped nuclear structure
/// for the heat transfer module
pub mod lumped_nuclear_structure_inspired_functions;
use gen_foam_lumped_nuclear_structure::*;

/// this module contains constructors for an array cv which 
/// functions as a one dimensional cartesian conduction medium
pub mod one_dimension_cartesian_conducting_medium;
pub use one_dimension_cartesian_conducting_medium::*;


use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::{HeatTransferInteractionType, enum_selection_alpha::single_control_vol_interactions::calculate_between_two_singular_cv_nodes};

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::
enum_selection_alpha::
single_control_volume_timestep_control::
calculate_mesh_stability_conduction_timestep_for_single_node_and_bc;

use crate::heat_transfer_lib::thermophysical_properties::Material::*;



use super::{ArrayCVType, SingleCVNode};



// Solve `Ax=b`
fn solve_example_array() -> Result<(), error::LinalgError> {
    let a: Array2<f64> = 
    arr2(&[[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]]);
    let b: Array1<f64> = arr1(&[1.0,2.0,3.01]);
    let _x = a.solve(&b)?;
    Ok(())
}

// Solve `Ax=b` for many b with fixed A
fn factorize() -> Result<(), error::LinalgError> {
    let a: Array2<f64> = random((3, 3));
    let f = a.factorize_into()?; // LU factorize A (A is consumed)
    for _ in 0..10 {
        let b: Array1<f64> = random(3);
        let _x = f.solve_into(b)?; // solve Ax=b using factorized L, U
    }
    Ok(())
}

// sandbox 
fn sandbox() -> Result<(), error::LinalgError> {

    // i want to index into the array and change stuff
    let mut a: Array2<f64> = 
    arr2(&[[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]]);


    // https://docs.rs/ndarray/latest/ndarray/doc/ndarray_for_numpy_users/index.html
    //
    // the syntax to access by element is slightly different from c
    a[[1,2]] = 0.5;



    let b: Array1<f64> = arr1(&[1.0,2.0,3.01]);
    let _x = a.solve(&b)?;


    // we can also initialise an array using 
    // 
    // surprisingly, we can use thermal conductance as well
    // and hopefully in general, 

    let mut thermal_conductance_matrix: Array2<ThermalConductance> = Array::zeros((3,3));

    thermal_conductance_matrix[[0,0]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    thermal_conductance_matrix[[0,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[0,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    thermal_conductance_matrix[[1,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[1,1]] = ThermalConductance::new::<watt_per_kelvin>(3.0);
    thermal_conductance_matrix[[1,2]] = ThermalConductance::new::<watt_per_kelvin>(1.0);

    thermal_conductance_matrix[[2,0]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[2,1]] = ThermalConductance::new::<watt_per_kelvin>(1.0);
    thermal_conductance_matrix[[2,2]] = ThermalConductance::new::<watt_per_kelvin>(3.0);

    // this is the power vector
    let mut power_vector: Array1<Power> = Array::zeros(3);

    power_vector[0] = Power::new::<watt>(1.0);
    power_vector[1] = Power::new::<watt>(2.0);
    power_vector[2] = Power::new::<watt>(3.01);


    // solve linear system
    let _temperature_vector: Array1<ThermodynamicTemperature> = 
    solve_conductance_matrix_power_vector(thermal_conductance_matrix,
        power_vector)?;


    // now how to solve this? I'll probably have to convert this into 
    // a float, but I lose my unit safety
    //
    // we can do type conversion here: 
    // https://github.com/rust-ndarray/ndarray/blob/master/examples/type_conversion.rs
    //
    // direct solving is problematic

    //let T = M.solve(&S)?;
    //
    // this is because ndarray allows you to create non dimensional 
    // arrays of anytime but ndarray-linalg only deals with f64
    // or other numeric types
    //
    // one other way is to convert all into f64, and verify the solution 
    // by mapping it back to its appropriate units 
    //
    // so T will be an f64 vector, which should be temperature 
    // S will be an f64 vector as well 
    // M will be an f64 matrix 
    //
    // so in effect, a custom solve method with appropriate unit 
    // checks would suffice
    //

    Ok(())
}

/// contains implementations for arrayCVType 
impl ArrayCVType {

    /// sets the current temperature vector and other 
    /// properties based on temperatures calculated
    /// for the next timestep
    /// and also other cleanup work
    pub fn advance_timestep(&mut self) -> Result<(),String>{
        return Err("not implemented".to_string());
    }
    /// gets the bulk temperature for the ArrayCV
    pub fn get_bulk_temperature(&mut self) -> 
    Result<ThermodynamicTemperature,String>{
        return Err("not implemented".to_string());
    }

    /// gets the maximum timestep for the arrayCV
    pub fn get_max_timestep(&mut self) -> Result<Time, String>{
        return Err("not implemented".to_string());
    }

    /// attaches a single cv to the front,entrance,
    /// lower or inner side of the 
    /// array cv
    /// 
    ///
    /// basically in whatever coordinate system, it is the lowest 
    /// value 
    ///
    /// for spheres, the lowest r (inner side)
    /// for cylinders, the lowest r or z (inner or lower)
    /// for cartesian, the lowest x, y or z (back)
    ///
    /// We use this convention because heat can flow either way 
    /// and fluid can flow either way. 
    ///
    /// So it's better to define higher and lower based upon a coordinate 
    /// axis
    pub fn link_single_cv_to_lower_side(&mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<(), String>{


        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.inner_single_cv
            },
        };

        // now link both cvs or calculate between them

        calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches a single cv to the exit,back,
    /// higher or outer side of the 
    /// array cv
    /// 
    ///
    /// basically in whatever coordinate system, it is the lowest 
    /// value 
    ///
    /// for spheres, the highest r (outer side)
    /// for cylinders, the highest r or z (outer or higher)
    /// for cartesian, the highest x, y or z (front)
    ///
    /// We use this convention because heat can flow either way 
    /// and fluid can flow either way. 
    ///
    /// So it's better to define higher and lower based upon a coordinate 
    /// axis
    pub fn link_single_cv_to_higher_side(&mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<(), String>{

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.outer_single_cv
            },
        };

        // now link both cvs or calculate between them

        calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches an array control volume to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    pub fn link_array_cv_to_the_front_of_this_cv(
        &mut self,
        array_cv_other: &mut ArrayCVType,
        interaction: HeatTransferInteractionType,) -> Result<(), String>{

        // basically we need to get the front of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.outer_single_cv
            },
        };

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        match array_cv_other {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.inner_single_cv
            },
        };

        calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches an constant heat flux BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant heat flux)
    pub fn link_heat_flux_bc_to_front_of_this_cv(
        &mut self,
        heat_flux_into_control_vol: HeatFluxDensity,
        interaction: HeatTransferInteractionType) -> Result<(),String>{

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
        };

        // once area is calculated, we can calculate heat flowrate into 
        // cv
        let heat_flowrate_into_control_vol: Power = 
        heat_flux_into_control_vol * heat_transfer_area;

        // then push the power to the front control volume, 

        match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                cartesian_array_cv.
                    outer_single_cv.rate_enthalpy_change_vector 
                    .push(heat_flowrate_into_control_vol);

                // needed to check the timescales between the 
                // outer or front surface node to the bc
                let cv_material = cartesian_array_cv.outer_single_cv.material_control_volume;
                match cv_material {
                    Solid(_) => {
                        // in this case, we just have one cv and one bc 
                        // so we only consider thermal inertia of this cv 
                        calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                            &mut cartesian_array_cv.outer_single_cv,
                            interaction)?;
                        ()
                    },
                    Liquid(_) => {
                        todo!("need to calculate convection based time scales")
                    },
                }
            },
        }

        return Ok(());
    }

    /// attaches an constant heat flux BC to the front of this 
    /// array control volume 
    /// (constant heat flux) ---- (back --- cv_self --- front)
    pub fn link_heat_flux_bc_to_back_of_this_cv(
        &mut self,
        heat_flux_into_control_vol: HeatFluxDensity,
        interaction: HeatTransferInteractionType) -> Result<(),String>{

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
        };

        // once area is calculated, we can calculate heat flowrate into 
        // cv
        let heat_flowrate_into_control_vol: Power = 
        heat_flux_into_control_vol * heat_transfer_area;

        // then push the power to the front control volume, 

        match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                cartesian_array_cv.
                    inner_single_cv.rate_enthalpy_change_vector 
                    .push(heat_flowrate_into_control_vol);

                // needed to check the timescales between the 
                // inner or front surface node to the bc
                let cv_material = cartesian_array_cv.inner_single_cv.material_control_volume;
                match cv_material {
                    Solid(_) => {
                        // in this case, we just have one cv and one bc 
                        // so we only consider thermal inertia of this cv 
                        calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
                            &mut cartesian_array_cv.inner_single_cv,
                            interaction)?;
                        ()
                    },
                    Liquid(_) => {
                        todo!("need to calculate convection based time scales")
                    },
                }
            },
        }

        return Ok(());
    }

    /// attaches an constant heat rate BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant heat rate)
    pub fn link_heat_addition_to_front_of_this_cv(
        &mut self,
        heat_rate: Power,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        Err("array cv not yet implemented".to_string())
    }

    /// attaches an constant heat rate BC to the front of this 
    /// array control volume 
    /// (constant heat rate) --- (back --- cv_self --- front)
    pub fn link_heat_addition_to_back_of_this_cv(
        &mut self,
        heat_rate: Power,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        Err("array cv not yet implemented".to_string())
    }

    /// attaches an constant temperature BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant temperature)
    pub fn link_constant_temperature_to_front_of_this_cv(
        &mut self,
        bc_temperature: ThermodynamicTemperature,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        Err("array cv not yet implemented".to_string())
    }

    /// attaches an constant temperature BC to the front of this 
    /// array control volume 
    /// (constant temperature) --- (back --- cv_self --- front)
    pub fn link_constant_temperature_to_back_of_this_cv(
        &mut self,
        bc_temperature: ThermodynamicTemperature,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        Err("array cv not yet implemented".to_string())
    }

    /// calculates timestep for a single cv attached to the front of the 
    /// array cv
    /// (back --- cv_self --- front) ---- (single cv)
    pub fn calculate_timestep_for_single_cv_to_front_of_array_cv(
        &mut self,
        single_cv_node: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        Err("array cv not yet implemented".to_string())
    }

    /// calculates timestep for a single cv attached to the back of the 
    /// array cv
    /// (single cv) --- (back --- cv_self --- front)
    pub fn calculate_timestep_for_single_cv_to_back_of_array_cv(
        &mut self,
        single_cv_node: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        Err("array cv not yet implemented".to_string())
    }


    /// calculates timestep for an array cv attached to the front of the 
    /// array cv
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    pub fn calculate_timestep_for_array_cv_to_front_of_this_array_cv(
        &mut self,
        array_cv_other: &mut ArrayCVType,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        Err("array cv not yet implemented".to_string())
    }
}
