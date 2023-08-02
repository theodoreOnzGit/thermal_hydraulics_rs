use std::f64::consts::PI;

use peroxide::prelude::P;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::CVType::*;


use crate::heat_transfer_lib::control_volume_calculations
::heat_transfer_interactions::*;

use crate::heat_transfer_lib::thermophysical_properties::Material
::{Solid,Liquid};

use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_interactions::enum_selection_alpha::*;
use crate::heat_transfer_lib::thermophysical_properties::thermal_diffusivity::thermal_diffusivity;

/// calculates an appropriate time step for two single control volume 
/// nodes, appending the max time step allowable for each one
/// first based on thermal conductivity
///
/// then courant number (TBC)
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_mesh_stability_timestep_for_two_single_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) 
-> Result<ThermalConductance, String> 
{

    let temperature_1: ThermodynamicTemperature = 
    single_cv_1.get_temperature()?;

    let temperature_2: ThermodynamicTemperature = 
    single_cv_2.get_temperature()?;

    let pressure_1: Pressure = 
    single_cv_1.pressure_control_volume.clone();

    let pressure_2: Pressure = 
    single_cv_2.pressure_control_volume.clone();

    let conductance: ThermalConductance = match 
        interaction {
            HeatTransferInteractionType::UserSpecifiedThermalConductance(
                user_specified_conductance) => user_specified_conductance,
            HeatTransferInteractionType
                ::SingleCartesianThermalConductanceOneDimension(
                material,thickness) => get_conductance_single_cartesian_one_dimension(
                    material,
                    temperature_1, 
                    temperature_2, 
                    pressure_1, 
                    pressure_2, 
                    thickness)?,
            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(
                (solid_material, shell_thickness,
                solid_temperature, solid_pressure), 
                (h, inner_diameter, cylinder_length)) => {

                    let id: Length = inner_diameter.clone().into();
                    let thicnkess: Length = shell_thickness.clone().into();

                    let od: Length = id+thicnkess;

                    let outer_diameter: OuterDiameterThermalConduction = 
                    OuterDiameterThermalConduction::from(od);

                    // after all the typing conversion, we can 
                    // get our conductance
                    get_conductance_single_cylindrical_radial_solid_liquid(
                        solid_material,
                        solid_temperature,
                        solid_pressure,
                        h,
                        inner_diameter,
                        outer_diameter,
                        cylinder_length,
                        CylindricalAndSphericalSolidFluidArrangement::
                        FluidOnInnerSurfaceOfSolidShell ,
                    )?
                },

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidOutside(
                (solid_material, shell_thickness,
                solid_temperature, solid_pressure), 
                (h, outer_diameter, cylinder_length)) => {

                    let od: Length = outer_diameter.clone().into();
                    let thicnkess: Length = shell_thickness.clone().into();

                    let id: Length = od - thicnkess;

                    let inner_diameter: InnerDiameterThermalConduction = 
                    InnerDiameterThermalConduction::from(id);

                    // after all the typing conversion, we can 
                    // get our conductance
                    get_conductance_single_cylindrical_radial_solid_liquid(
                        solid_material,
                        solid_temperature,
                        solid_pressure,
                        h,
                        inner_diameter,
                        outer_diameter,
                        cylinder_length,
                        CylindricalAndSphericalSolidFluidArrangement::
                        FluidOnInnerSurfaceOfSolidShell ,
                    )?
                },
            // note: actually function signatures are a little more 
            // friendly to use than packing enums with lots of stuff 
            // so may change stuffing enums with tuples to stuffing 
            // enums with a single struct
            HeatTransferInteractionType
                ::DualCylindricalThermalConductance(
                (inner_material,inner_shell_thickness),
                (outer_material,outer_shell_thickness),
                (inner_diameter,
                outer_diameter,
                cylinder_length)
            ) => {
                    // first, want to check if inner_diameter + 
                    // shell thicknesses is outer diameter 

                    let expected_outer_diameter: Length;
                    let id: Length = inner_diameter.into();
                    let inner_thickness: Length =  inner_shell_thickness.into();
                    let outer_thickness: Length =  outer_shell_thickness.into();

                    expected_outer_diameter = 
                        id + inner_thickness + outer_thickness;

                    let od: Length = outer_diameter.into();

                    // inner diameter and outer diameter values must be 
                    // equal to within 1 nanometer 1e-9 m
                    if (od.value - expected_outer_diameter.value).abs() > 1e-9
                    {

                        let mut error_str: String = "the inner diameter 
                            plus shell thicknesses do not equate 
                            to outer diameter".to_string();

                        error_str += "supplied outer diameter (m):";
                        error_str += &od.value.to_string();
                        error_str += "expected outer diameter (m):";
                        error_str += &expected_outer_diameter.value.to_string();


                        return Err(error_str
                        );
                    }

                    get_conductance_cylindrical_radial_two_materials(
                        inner_material,
                        outer_material,
                        temperature_1, //convention, 1 is inner shell
                        temperature_2, // convention 2, is outer shell
                        pressure_1,
                        pressure_2,
                        inner_diameter,
                        inner_shell_thickness,
                        outer_shell_thickness,
                        cylinder_length,
                    )?
                },
            HeatTransferInteractionType::UserSpecifiedHeatAddition  
                => {
                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());
                },

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCustomArea(_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(_,_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

                },
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalInnerArea(_,_) => {

                    return Err("interaction type needs to be \n 
                        thermal conductance".to_string());

                },
            HeatTransferInteractionType::
                DualCartesianThermalConductance(
                (material_1, thickness_1),
                (material_2,thickness_2)) => { 
                    
                    let conductnace_layer_1: ThermalConductance 
                    = get_conductance_single_cartesian_one_dimension(
                        material_1,
                        temperature_1, 
                        temperature_1, 
                        pressure_1, 
                        pressure_1, 
                        thickness_1)?;

                    let conductnace_layer_2: ThermalConductance 
                    = get_conductance_single_cartesian_one_dimension(
                        material_2,
                        temperature_2, 
                        temperature_2, 
                        pressure_2, 
                        pressure_2, 
                        thickness_2)?;

                    let overall_resistance = 
                    1.0/conductnace_layer_2 
                    + 1.0/conductnace_layer_1;

                    // return the conductance or resistnace inverse

                    1.0/overall_resistance
            },
            HeatTransferInteractionType::
                DualCartesianThermalConductanceThreeDimension(
                data_dual_cartesian_conduction) 
                => {

                    let material_1 = 
                    data_dual_cartesian_conduction .material_1;

                    let material_2 = 
                    data_dual_cartesian_conduction .material_2;

                    let thickness_1 = 
                    data_dual_cartesian_conduction .thickness_1;

                    let thickness_2 = 
                    data_dual_cartesian_conduction .thickness_2;

                    let xs_area = 
                    data_dual_cartesian_conduction .xs_area;

                    get_conductance_dual_cartesian_three_dimensions(
                        material_1, 
                        material_2, 
                        temperature_1, 
                        temperature_2, 
                        pressure_1, 
                        pressure_2, 
                        xs_area, 
                        thickness_1,
                        thickness_2)?
                },

            HeatTransferInteractionType::
                UserSpecifiedConvectionResistance(
                data_convection_resistance) 
                => {

                    let heat_transfer_coeff: HeatTransfer = 
                    data_convection_resistance.heat_transfer_coeff;
                    let surf_area: Area = 
                    data_convection_resistance.surf_area.into();

                    heat_transfer_coeff * surf_area
                },
        };

    return Ok(conductance);
}


/// Calculates a time step based on fourier number for a single 
/// control volume, usually tied to some boundary condition
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
    control_vol: &mut SingleCVNode,
    interaction: HeatTransferInteractionType) -> Result<Time,String> {

    // here we have timestep based on the generic lengthscale of the 
    // control volume 
    let mut cv_timestep:Time = 
    control_vol.calculate_conduction_timestep()?;

    // we may have other time scales based on differing length scales 
    // of the control volume 
    //
    // so we will need to calculate time scales based on these other 
    // length scales and then calculate each of their own time scales.
    // if shorter, then we need to append it to the respective control 
    // volumes

    let cv_material = control_vol.material_control_volume.clone();
    let cv_pressure = control_vol.pressure_control_volume.clone();
    let cv_temperature = control_vol.get_temperature()?;

    
    let cv_alpha: DiffusionCoefficient = 
    thermal_diffusivity(cv_material,
        cv_temperature,
        cv_pressure)?;


    let max_mesh_fourier_number: f64 = 0.25;


    match interaction {
        HeatTransferInteractionType::
            UserSpecifiedHeatAddition => {

                // do nothing

                ()
            },
        HeatTransferInteractionType::
            UserSpecifiedThermalConductance(_) => {

                // if a conductance is specified, don't 
                // do anything

            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(area) => {
                // when a normal area is given,
                // we can calculate volume to area ratio 

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }

            },

        HeatTransferInteractionType::
            SingleCartesianThermalConductanceOneDimension(
            material, x_thickness) => {

                // the given material here overrides the normal 
                // material 
                let cv_alpha: DiffusionCoefficient = 
                thermal_diffusivity(material,
                    cv_temperature,
                    cv_pressure)?;

                // if we connect the cv to a boundary condition,
                // then the length provided here is what we need 
                // to bother with 

                let lengthscale: Length = x_thickness.into();

                let time_step_max_based_on_x_thickness: Time 
                = max_mesh_fourier_number *
                lengthscale * 
                lengthscale / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_x_thickness {
                    cv_timestep = time_step_max_based_on_x_thickness;
                }
            },

        HeatTransferInteractionType::
            DualCartesianThermalConductance(
            (material_1, length_1), 
            (material_2,length_2)) => {
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;

                let length_1: Length = length_1.into();
                let length_2: Length = length_2.into();

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                length_2 *
                length_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()
            },

        HeatTransferInteractionType::
            DualCylindricalThermalConductance(
            (material_1, radius_1), 
            (material_2,radius_2), 
            _) => {
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more radiuss
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;

                let radius_1: Length = radius_1.into();
                let radius_2: Length = radius_2.into();

                let timestep_1: Time = max_mesh_fourier_number * 
                radius_1 *
                radius_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                radius_2 *
                radius_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()
            },


        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(
            data_dual_cartesian_conduction_data) => {

                let material_1 = data_dual_cartesian_conduction_data.
                    material_1.clone();

                let material_2 = data_dual_cartesian_conduction_data.
                    material_2.clone();


                let length_1 : Length = data_dual_cartesian_conduction_data.
                    thickness_1.clone().into();

                let length_2 : Length = data_dual_cartesian_conduction_data.
                    thickness_2.clone().into();
                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more radiuss
                //
                // you only have one control volume bascially,
                // and you should only use dual cartesian thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep
                //
                // the other consideration is to take the shorter of 
                // the two time steps and put it into the cv timestep 
                let alpha_1: DiffusionCoefficient = 
                thermal_diffusivity(material_1,
                    cv_temperature,
                    cv_pressure)?;

                let alpha_2: DiffusionCoefficient = 
                thermal_diffusivity(material_2,
                    cv_temperature,
                    cv_pressure)?;


                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                let timestep_2: Time = max_mesh_fourier_number * 
                length_2 *
                length_2 / 
                alpha_2;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }

                if cv_timestep > timestep_2 {
                    cv_timestep = timestep_2;
                }

                // done!
                ()

            },

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
            (material,radius,
            temperature,pressure),_) => {

                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cylindrical thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep

                let alpha_1: DiffusionCoefficient = 
                thermal_diffusivity(material,
                    temperature,
                    pressure)?;

                let length_1: Length =  radius.into();

                let length_1 = length_1*0.5;

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }
                ()

            },

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(
            (material,radius,
            temperature,pressure),_) => {

                // for a single node connected to a BC, you're 
                // not really supposed to have a timescale 
                // based on two or more lengths
                //
                // you only have one control volume bascially,
                // and you should only use dual cylindrical thermal 
                // conductance for two control volumes
                // I won't do anything based on this 
                // or just use the generic timestep

                let alpha_1: DiffusionCoefficient = 
                thermal_diffusivity(material,
                    temperature,
                    pressure)?;

                let length_1: Length =  radius.into();

                let length_1 = length_1*0.5;

                let timestep_1: Time = max_mesh_fourier_number * 
                length_1 *
                length_1 / 
                alpha_1;

                if cv_timestep > timestep_1 {
                    cv_timestep = timestep_1;
                }
                ()

            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(l, od) => {

                // this is treated like a custom area kind of thing 
                // so we calculate the area first

                let cylinder_length: Length = l.into();
                let outer_diameter: Length = od.into();

                let area: Area = PI * outer_diameter * cylinder_length;

                // and then do the boilerplate code

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(l, id) => {

                // this is treated like a custom area kind of thing 
                // so we calculate the area first

                let cylinder_length: Length = l.into();
                let inner_diameter: Length = id.into();

                let area: Area = PI * inner_diameter * cylinder_length;

                // and then do the boilerplate code

                let cv_volume = control_vol.volume.clone();

                let volume_to_area_ratio: Length = cv_volume/area;

                // we can calculate a timestep

                let time_step_max_based_on_volume_to_area: Time 
                = max_mesh_fourier_number *
                volume_to_area_ratio * 
                volume_to_area_ratio / 
                cv_alpha;

                // if the max timestep is shorter than this calculated 
                // cv timestep, use it

                if cv_timestep > time_step_max_based_on_volume_to_area {
                    cv_timestep = time_step_max_based_on_volume_to_area;
                }
            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(_) => {

                // if a resistance is specified, don't 
                // do anything

                ()
            },

    }


    control_vol.max_timestep_vector.push(cv_timestep);

    return Ok(cv_timestep);

}
