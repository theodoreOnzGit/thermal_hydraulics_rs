use std::f64::consts::PI;

use uom::si::thermal_conductance::watt_per_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib::control_volume_calculations::common_functions::obtain_thermal_conductance_annular_cylinder;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use crate::heat_transfer_lib::thermophysical_properties ::Material;

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::calculations::
CylindricalAndSphericalSolidFluidArrangement::*;

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
/// 
/// Now, it is important also to specify which control volume is 
/// adjacent to the
/// the inner radius and which one is at the outer radius
/// 
///
pub(in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions) 
fn get_conductance_cylindrical_radial_two_materials(
    material_inner_shell: Material,
    material_outer_shell: Material,
    material_temperature_inner_shell: ThermodynamicTemperature,
    material_temperature_outer_shell: ThermodynamicTemperature,
    material_pressure_inner_shell: Pressure,
    material_pressure_outer_shell: Pressure,
    id: InnerDiameterThermalConduction,
    inner_shell_thickness: RadialCylindricalThicknessThermalConduction,
    outer_shell_thickness: RadialCylindricalThicknessThermalConduction,
    l: CylinderLengthThermalConduction) -> Result<ThermalConductance,String> 
{

    // convert quantities to their standard uom Length types
    let id: Length = id.into();
    let inner_shell_thickness: Length = inner_shell_thickness.into();
    let outer_shell_thickness: Length = outer_shell_thickness.into();
    let od: Length = id + inner_shell_thickness + outer_shell_thickness;
    let l: Length = l.into();


    // also, we need the interim_diameter which is the outer 
    // diameter of the inner shell 
    // so that we know the relative thicknesses of these two layers
    let interim_diameter: Length = id + inner_shell_thickness;

    // now we need to calculate thermal resistance of an annular 
    // cylindrical layer

    // so we have inner layer thermal conductance 
    // the question mark propagates the error with match statement

    let inner_thermal_conductivity = thermal_conductivity(
        material_inner_shell,
        material_temperature_inner_shell,
        material_pressure_inner_shell)?;

    let outer_thermal_conductivity = thermal_conductivity(
        material_outer_shell,
        material_temperature_outer_shell,
        material_pressure_outer_shell)?;

    // again 
    // |                                 | 
    // |                                 | 
    // *---------------*-----------------* 
    // |                                 | 
    // |                                 | 
    //
    // inner        interim             outer 
    //
    let inner_layer_conductance: ThermalConductance = 
    obtain_thermal_conductance_annular_cylinder(
        id,
        interim_diameter,
        l,
        inner_thermal_conductivity)?;

    let outer_layer_conductance: ThermalConductance = 
    obtain_thermal_conductance_annular_cylinder(
        interim_diameter,
        od,
        l,
        outer_thermal_conductivity)?;

    // the obtain_thermal_conductance_annular_cylinder already checks 
    // that conductivity should be non zero, 
    // k should be non zero and lengths are non zero
    // also checks that 
    // od > id

    // now we have both conductances, we can sum up their resistances 

    let total_thermal_resistance = 
    1.0/inner_layer_conductance + 1.0/outer_layer_conductance;

    return Ok(1.0/total_thermal_resistance);


}



/// Suppose we have two control volumes of differing materials and  
/// temperature 
///
/// one control volume is a solid and the other is a fluid
///
/// now, we want to calculate thermal resistance between them
///
/// for fluids, the thermal resistance or conductance is quite 
/// straightforward 
///
/// from fluid to solid heat transfer,
/// Q = -hA (T_solid - T_fluid)
///
/// resistance here is hA 
/// where A is the curved surface area pi*D*L
///
/// for solid thermal resistance, we use the 
/// obtain_thermal_conductance_annular_cylinder 
/// function under common functions
///
/// that would need an inner diameter and an outer diameter
///
/// There are two cases here. 
/// 
/// Firstly, 
/// the fluid is an in the tube side of a heat exchanger or pipe,
/// hence the solid is considered on the outside 
///
/// The surface area will be based on the inner diameter 
///
/// Secondly, the fluid is on the outside of the cylindrical solid,
/// in this case, the surface area will be based on the outer diameter
///
/// 
///
///
pub(in crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions) 
fn get_conductance_single_cylindrical_radial_solid_liquid(
    solid: Material,
    solid_temperature: ThermodynamicTemperature,
    solid_pressure: Pressure,
    h: HeatTransfer,
    id: InnerDiameterThermalConduction,
    od: OuterDiameterThermalConduction,
    l: CylinderLengthThermalConduction,
    solid_liquid_arrangement: CylindricalAndSphericalSolidFluidArrangement) 
-> Result<ThermalConductance,String> 
{

    // first thing first, let's do the uncomplicated part 
    // that is to calculate the thermal resistance of the solid part

    // convert quantities to their standard uom Length types
    let id: Length = id.into();
    let od: Length = od.into();
    let l: Length = l.into();

    // I need to typecast the solid material into a material enum 


    let solid_thermal_conductivity: ThermalConductivity = 
    thermal_conductivity(solid, 
        solid_temperature, 
        solid_pressure)?;

    let solid_layer_conductance: ThermalConductance = 
    obtain_thermal_conductance_annular_cylinder(
        id,
        od,
        l,
        solid_thermal_conductivity)?;

    // now then, we deal with the thermal resistance on the liquid 
    // end
    //
    // firstly, we get a h value, 
    // which is already in the function inputs
    // next, we need something to tell us if the solid is on the 
    // inside or outside, I could either use an enum or an if statement
    //
    // enums are strongly typed so I'll use those ,
    // this way, the compiler will guide the user on what to do
    //
    // it's kind of long but the compiler will tell you what inputs 
    // you need basically

    // now let's handle the cases, and depending on the cases, 
    // the surface area will change 
    //
    // area is PI * D * L


    let surface_area_for_solid_liquid_boundary: Area = match 
        solid_liquid_arrangement {
            FluidOnInnerSurfaceOfSolidShell => PI * id * l,
            FluidOnOuterSurfaceOfSolidShell => PI * od * l,
        };

    let liquid_thermal_conductance: ThermalConductance = 
    h * surface_area_for_solid_liquid_boundary;

    // now we do overall thermal resistance and conductance

    let overall_thermal_resistance = 1.0/liquid_thermal_conductance 
    + 1.0/solid_layer_conductance;
    
    return Ok(1.0/overall_thermal_resistance);

}

/// basically an enum for you to specify 
/// if the liquid on the inner curved surface of the shell or outer 
/// curved surface of the shell
///
/// in the context of a convection and conductivity 
/// thermal resistance calculation,
///
pub enum CylindricalAndSphericalSolidFluidArrangement {
    FluidOnInnerSurfaceOfSolidShell,
    FluidOnOuterSurfaceOfSolidShell
}


