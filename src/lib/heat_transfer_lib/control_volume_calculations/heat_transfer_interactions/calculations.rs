use std::f64::consts::PI;

use uom::si::thermal_conductance::watt_per_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{f64::*, pressure::atmosphere};
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use crate::heat_transfer_lib::thermophysical_properties
::{Material, specific_enthalpy::specific_enthalpy};
use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::HeatTransferEntity;

/// for 1D calculations, we need to calculate conductance as well,
/// but there is no area, hence, we have to use a unit area to calculate 
/// the conductance 
pub const UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS: f64 = 1.0;

/// Suppose we have two control volumes of the same materials and  
/// temperature and we put a 1D thermal resistance between them 
///
/// we would need to return a thermal conductance based on a 1D 
/// heat transfer model
///
/// conductance is watts per kelvin or 
/// q = (kA)/dx * dT
/// conductance here is kA/dx
/// thermal resistance is 1/conductance
///
/// For a 1D case, the area is not defined, but I'm giving it a unit 
/// area value of 1 meter squared specific to 1D calculations
///
/// Note that the control volume MUST have the same cross sectional 
/// area so that it is consistent. Will need to be a 1D control 
/// volume of sorts
///
pub(in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions) 
fn get_conductance_single_cartesian_one_dimension(
    material: Material,
    material_temperature_1: ThermodynamicTemperature,
    material_temperature_2: ThermodynamicTemperature,
    material_pressure_1: Pressure,
    material_pressure_2: Pressure,
    thickness: XThicknessThermalConduction) -> Result<ThermalConductance,String> 
{

    // note, the question mark here (?) at the end denotes 
    // error propagation.
    //
    // So basically, we do pattern matching. If it's an error,
    // return an error as the result enum 
    // if it's successful, it will obtain the value in the Ok()
    // variant of the enum

    let average_material_pressure: Pressure = 
    (material_pressure_1 + material_pressure_2)/2.0;

    let average_material_temperature_kelvin = 
    material_temperature_1.value + material_temperature_2.value
    -273.15;

    let average_material_temperature = 
    ThermodynamicTemperature::new::<kelvin>(average_material_temperature_kelvin);


    let material_thermal_conductivity = thermal_conductivity(
        material, 
        average_material_temperature, 
        average_material_pressure)?;    

    let cross_sectional_area: Area = Area::new::<uom::si::area::square_meter>(
        UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS);

    let material_thickness: Length = thickness.into();

    let material_conductance: ThermalConductance = 
    material_thermal_conductivity * 
    cross_sectional_area / 
    material_thickness;
    

    return Ok(material_conductance);

}

/// Suppose we have two control volumes of differing materials and  
/// temperature and we put a thermal resistance between them 
/// in the x coordinate
///
/// we would need to return a thermal conductance based on a 1D 
/// heat transfer model
///
/// conductance is watts per kelvin or 
/// q = (kA)/dx * dT
/// conductance here is kA/dx
/// thermal resistance is 1/conductance
/// 
///
///
pub(in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions) 
fn get_conductance_single_cartesian_three_dimensions(
    material_1: Material,
    material_2: Material,
    material_temperature_1: ThermodynamicTemperature,
    material_temperature_2: ThermodynamicTemperature,
    material_pressure_1: Pressure,
    material_pressure_2: Pressure,
    xs_area: CrossSectionalArea,
    thickness: XThicknessThermalConduction) -> Result<ThermalConductance,String> 
{

    // note, the question mark here (?) at the end denotes 
    // error propagation.
    //
    // So basically, we do pattern matching. If it's an error,
    // return an error as the result enum 
    // if it's successful, it will obtain the value in the Ok()
    // variant of the enum



    let material_1_thermal_conductivity = thermal_conductivity(
        material_1, 
        material_temperature_1, 
        material_pressure_1)?;    


    let material_2_thermal_conductivity = thermal_conductivity(
        material_2, 
        material_temperature_2, 
        material_pressure_2)?;    

    let material_thickness: Length = thickness.into();

    let cross_sectional_area: Area = xs_area.into();

    // we then calculate conductances of each material

    let material_1_conductance: ThermalConductance = 
    material_1_thermal_conductivity * 
    cross_sectional_area / 
    material_thickness;

    let material_2_conductance: ThermalConductance = 
    material_2_thermal_conductivity * 
    cross_sectional_area / 
    material_thickness;

    // we invert them to find the resistance
    //
    // However, doing inversion is problematic if either conductance 
    // is zero (infinite resistance)
    //
    // so if either conductance is zero, then we just return a zero 
    // conductance

    if material_1_conductance.value == 0.0 || 
        material_2_conductance.value == 0.0 {
        return Ok(ThermalConductance::new::<watt_per_kelvin>(0.0));
    }

    // if less than zero, then we have an issue

    if material_1_conductance.value < 0.0 || 
        material_2_conductance.value < 0.0 {
        return Err("thermal negative thermal conductance".to_string());
    }


    let overall_thermal_resistance = 
    1.0/material_2_conductance 
    + 1.0/material_1_conductance;

    let overall_conductance: ThermalConductance = 
    1.0/overall_thermal_resistance;
    

    return Ok(overall_conductance);

}

/// Suppose we have two control volumes of differing materials and  
/// temperature and we put a thermal resistance between them 
/// in the cylindrical radial coordinate
///
/// we would need to return a thermal conductance based on a 1D 
/// heat transfer model in the r coordinate
///
/// conductance is watts per kelvin or 
/// q = (kA)/dr * dT
/// conductance here is kA/dr
/// thermal resistance is 1/conductance
/// 
/// Now, it is important also to specify which control volume is 
/// adjacent to the
/// the inner radius and which one is at the outer radius
/// 
///
pub(in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions) 
fn get_conductance_single_cylindrical_radial(
    material_inner_shell: Material,
    material_outer_shell: Material,
    material_temperature_inner_shell: ThermodynamicTemperature,
    material_temperature_outer_shell: ThermodynamicTemperature,
    material_pressure_inner_shell: Pressure,
    material_pressure_outer_shell: Pressure,
    id: InnerDiameterThermalConduction,
    od: OuterDiameterThermalConduction,
    l: CylinderLengthThermalConduction,
    thickness: XThicknessThermalConduction) -> Result<ThermalConductance,String> 
{

    // convert quantities to their standard uom types
    let id: Length = id.into();
    let od: Length = od.into();
    let l: Length = l.into();

    // we can begin calculating surface areas
    let inner_surface_area: Area = PI * id * l;
    let outer_surface_area: Area = PI * od * l;

    todo!()


    //let material_thermal_conductivity = thermal_conductivity(
    //    material, 
    //    average_material_temperature, 
    //    average_material_pressure)?;    


}
