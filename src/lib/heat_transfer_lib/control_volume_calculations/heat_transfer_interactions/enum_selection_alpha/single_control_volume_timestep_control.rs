use std::f64::consts::PI;

use peroxide::prelude::P;
use uom::si::length::meter;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::f64::*;
use crate::heat_transfer_lib;
use crate::heat_transfer_lib::thermophysical_properties::LiquidMaterial;
use crate::heat_transfer_lib::thermophysical_properties::density::density;
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
use crate::heat_transfer_lib::thermophysical_properties::specific_heat_capacity::specific_heat_capacity;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
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
-> Result<Time, String> 
{

    let temperature_1: ThermodynamicTemperature = 
    single_cv_1.get_temperature()?;

    let temperature_2: ThermodynamicTemperature = 
    single_cv_2.get_temperature()?;

    let pressure_1: Pressure = 
    single_cv_1.pressure_control_volume.clone();

    let pressure_2: Pressure = 
    single_cv_2.pressure_control_volume.clone();

    let material_1 = single_cv_1.material_control_volume.clone();
    let material_2 = single_cv_2.material_control_volume.clone();

    // we use this to get the diffusion coefficient for both control 
    // volumes

    let alpha_1: DiffusionCoefficient = 
    thermal_diffusivity(material_1,
        temperature_1,
        pressure_1)?;

    let alpha_2: DiffusionCoefficient = 
    thermal_diffusivity(material_2,
        temperature_2,
        pressure_2)?;

    // we also want to get the ratios of the two thermal inertias 

    let cp_1: SpecificHeatCapacity = specific_heat_capacity(
        material_1,
        temperature_1,
        pressure_1)?;
    let cp_2: SpecificHeatCapacity = specific_heat_capacity(
        material_2,
        temperature_2,
        pressure_2)?;

    let mass_1: Mass = single_cv_1.mass_control_volume.clone();
    let mass_2: Mass = single_cv_2.mass_control_volume.clone();

    // we get heat capacity 

    let mcp_1: HeatCapacity = mass_1 * cp_1;
    let mcp_2: HeatCapacity = mass_2 * cp_2;

    // and the ratios of the heat capacities: 

    let heat_capacity_ratio_cv_1 = mcp_1 / (mcp_1 + mcp_2);
    let heat_capacity_ratio_cv_1: f64 = heat_capacity_ratio_cv_1.value;

    let heat_capacity_ratio_cv_2 = mcp_1 / (mcp_1 + mcp_2);
    let heat_capacity_ratio_cv_2: f64 = heat_capacity_ratio_cv_2.value;

    // the heat capacity ratios will be used to split up the length 
    // scales as needed

    // the strategy here is to find the shorter of the two lengthscales 
    // calculate the time step for each control volume, 
    // and then find the shortest timestep, then load both control volumes 
    // with the shortest time step
    let max_mesh_fourier_number: f64 = 0.25;

    // we'll calculate two baseline time scales
    let mut cv_1_timestep:Time = 
    single_cv_1.calculate_conduction_timestep()?;

    let mut cv_2_timestep:Time = 
    single_cv_2.calculate_conduction_timestep()?;

    // this next step estimates the minimum timestep based 
    // on the interaction type, and then loads both control volumes 
    // with this minimum timestep

    match interaction {
        HeatTransferInteractionType::UserSpecifiedThermalConductance(
            user_specified_conductance) => {

                let lengthscale_stability_vec_1 = 
                single_cv_1.mesh_stability_lengthscale_vector.clone();
                let lengthscale_stability_vec_2 = 
                single_cv_2.mesh_stability_lengthscale_vector.clone();

                let mut min_lengthscale = Length::new::<meter>(0.0);

                for length in lengthscale_stability_vec_1.iter() {
                    if *length > min_lengthscale {
                        min_lengthscale = *length;
                    }
                }
                
                for length in lengthscale_stability_vec_2.iter() {
                    if *length > min_lengthscale {
                        min_lengthscale = *length;
                    }
                }

                // now that we've gotten both length scales and gotten 
                // the shortest one, we'll need to obtain a timescale 
                // from the conductance 
                //
                // Delta t = Fo * rho * cp * Delta x * resistance
                //
                // we'll convert conductance into total resistance first 
                // 

                let total_resistance = 1.0/user_specified_conductance;

                // then we'll approximate the thermal resistance of 
                // cv 1 using the ratio of heat capacities 
                //
                // This is approximation, not exact, I'll deal with 
                // the approximation as I need to later.

                let approx_thermal_resistance_1= total_resistance *
                heat_capacity_ratio_cv_1;

                let approx_thermal_cond_1: ThermalConductance = 
                1.0/approx_thermal_resistance_1;

                // do note that thermal conductance is in watts per kelvin 
                // not watts per m2/K
                // area is factored in
                //
                // qA = H (T_1 - T_2)
                //
                // if i want HbyA or HbyV, then divide by the control volume
                //
                // H = kA/L for thermal conductivity case
                //  
                // bummer, no surface area!
                // I suppose I'll just use the minimum length scales 
                // then for conservative-ness

                let approx_thermal_conductivity_1: ThermalConductivity = 
                approx_thermal_cond_1 / min_lengthscale;



                let approx_thermal_resistance_2 = total_resistance *
                heat_capacity_ratio_cv_2;

                let approx_thermal_cond_2: ThermalConductance = 
                1.0/approx_thermal_resistance_2;

                let approx_thermal_conductivity_2: ThermalConductivity = 
                approx_thermal_cond_2 / min_lengthscale;

                // let's get timescale 1 and 2 estimate 

                let density_1: MassDensity 
                = density(material_1, temperature_1, pressure_1)?;

                let density_2: MassDensity 
                = density(material_2, temperature_2, pressure_2)?;

                let rho_cp_1: VolumetricHeatCapacity = density_1 * cp_1;
                let rho_cp_2: VolumetricHeatCapacity = density_2 * cp_2;

                let approx_alpha_1: DiffusionCoefficient = 
                approx_thermal_conductivity_1/rho_cp_1;

                let approx_alpha_2: DiffusionCoefficient = 
                approx_thermal_conductivity_2/rho_cp_2;


                // Fo  = alpha * Delta t / Delta x / Delta x 
                //
                // Delta t = Fo * Delta x * Delta x / alpha 
                //

                let timescale_1: Time = max_mesh_fourier_number 
                * min_lengthscale 
                * min_lengthscale 
                / approx_alpha_1;

                let timescale_2: Time = max_mesh_fourier_number 
                * min_lengthscale 
                * min_lengthscale 
                / approx_alpha_2;

                let interaction_minimum_timescale: Time; 
                
                if timescale_1 > timescale_2 {
                    interaction_minimum_timescale = timescale_2;
                } else {
                    interaction_minimum_timescale = timescale_1;
                }


                // take minimum of the timescales and assign it to 
                // the cv 1 or cv 2 timescales
                //

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }
                
                
            },
        HeatTransferInteractionType
            ::SingleCartesianThermalConductanceOneDimension(
            material,thickness) => {

                // let's get the thickness and split it among cv 1 
                // and cv 2 

                let conduction_thickness: Length = thickness.into();

                let thickness_cv_1: Length = 
                heat_capacity_ratio_cv_1 * 
                conduction_thickness;

                let thickness_cv_2: Length = 
                heat_capacity_ratio_cv_2 * 
                conduction_thickness;

                // Fo  = alpha * Delta t / Delta x / Delta x 
                //
                // Delta t = Fo * Delta x * Delta x / alpha 
                //

                let timescale_1: Time = max_mesh_fourier_number 
                * thickness_cv_1 
                * thickness_cv_1 
                / alpha_1;

                let timescale_2: Time = max_mesh_fourier_number 
                * thickness_cv_2 
                * thickness_cv_2 
                / alpha_2;


                let interaction_minimum_timescale: Time; 
                
                if timescale_1 > timescale_2 {
                    interaction_minimum_timescale = timescale_2;
                } else {
                    interaction_minimum_timescale = timescale_1;
                }


                // take minimum of the timescales and assign it to 
                // the cv 1 or cv 2 timescales
                //

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }


            },
        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidInside(
            (_solid_material, shell_thickness,
            _solid_temperature, _solid_pressure), 
            (h, inner_diameter, _cylinder_length)) => {

                // now for this, we are going to have to factor in 
                // h to get our nusselt number
                // unfortunately, we don't have a value for liquid 
                // thermal conductivity
                // not unless we match try matching the liquid 
                // material of either cv first 
                //
                // I'll just set it to material 1 by default

                let mut liquid_thermal_conductivity: ThermalConductivity
                = thermal_conductivity(
                    material_1,
                    temperature_1,
                    pressure_1)?;

                // we assume liquid is material 1, solid is material 
                // 2 
                // if we are wrong, swop it around
                let mut alpha_liquid: DiffusionCoefficient = 
                thermal_diffusivity(
                    material_1,
                    temperature_1,
                    pressure_1)?;
                
                let mut alpha_solid: DiffusionCoefficient = 
                thermal_diffusivity(
                    material_2,
                    temperature_2,
                    pressure_2)?;

                match material_1 {
                    Solid(_) => {
                        // if material 1 is solid, then we should obtain 
                        // liquid thermal conductivity from material 2
                        liquid_thermal_conductivity = thermal_conductivity(
                            material_2,
                            temperature_2,
                            pressure_2
                        )?;

                        alpha_solid = thermal_diffusivity(
                                material_1,
                                temperature_1,
                                pressure_1)?;

                        alpha_liquid = thermal_diffusivity(
                                material_2,
                                temperature_2,
                                pressure_2)?;
                        ()
                    },
                    Liquid(_) => {
                        // if it's liquid material, match it and obtain 
                        // thermal conductivity 
                        // if the liquid is in material 1, do nothing
                        ()
                    },
                }



                // if both material 1 and 2 are solid at the same time 
                // or liquid at the same time, this is WRONG
                //
                // it's a runtime error, not quite ideal
                // obviously, want a compile time error
                // but wait till future to debug this

                match (material_1, material_2) {
                    (Solid(_), Solid(_)) => {
                        return Err("should have 1 fluid and 1 solid".to_string());
                    },
                    (Liquid(_), Liquid(_)) => {
                        return Err("should have 1 fluid and 1 solid".to_string());
                    },
                    _ => (),
                }

                // todo: probably need to check which lengthscale is 
                // right for nusselt number 

                let id: Length = inner_diameter.clone().into();
                let nusselt_number_diameter: Ratio = 
                h * id / liquid_thermal_conductivity;

                let nusselt_number_diameter: f64 = 
                nusselt_number_diameter.value;

                // normally time scale is 
                // Delta t = Fo * Delta x * Delta x / alpha 
                // 
                // if convection at boundary is taken into account:
                // Delta t = Fo * Delta x * Delta x / alpha Nu_(delta x)
                //
                //
                // for cylindrical control vol on inside, 
                // let the lengthscale be 
                // the radius

                let liquid_lengthscale: Length = id/2.0;

                let nusselt_number_radius = nusselt_number_diameter/2.0;

                let liquid_timescale: Time = max_mesh_fourier_number * 
                    liquid_lengthscale * 
                    liquid_lengthscale / 
                    alpha_liquid / 
                    nusselt_number_radius;

                // the solid lengthscale shall be half the thickness 
                // of the wall tubing 

                let solid_lengthscale: Length = shell_thickness.into();


                let solid_lengthscale: Length = id/2.0;

                let solid_timescale: Time = max_mesh_fourier_number * 
                    solid_lengthscale * 
                    solid_lengthscale / 
                    alpha_solid ;

                // we find the minimum timescale and append it 
                // to both control volumes

                let interaction_minimum_timescale: Time; 
                
                if solid_timescale > liquid_timescale {
                    interaction_minimum_timescale = liquid_timescale;
                } else {
                    interaction_minimum_timescale = solid_timescale;
                }


                // take minimum of the timescales and assign it to 
                // the cv 1 or cv 2 timescales
                //

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }

            },

        HeatTransferInteractionType::
            CylindricalConductionConvectionLiquidOutside(
            (_solid_material, shell_thickness,
            _solid_temperature, _solid_pressure), 
            (h, _outer_diameter, _cylinder_length)) => {
                // 
                //
                // now for liquid on the outside, we really have no idea 
                // how much thermal inertia it has based on the 
                // lengthscale, so we don't have 
                // a clue as to its length scale directly
                // 
                // however, we can guestimate the lengthscale of the fluid 
                // volume by taking the lengthscale of the solid, and then 
                // scaling it by the ratio of the thermal inertia of 
                // the solid to the liquid 

                let thickness: Length = shell_thickness.clone().into();

                let solid_lengthscale = thickness;

                // to get liquid lengthscale, we must first find out 
                // which control volume contains the liquid 

                //
                // I'll just set it to material 1 by default

                let mut liquid_thermal_conductivity: ThermalConductivity
                = thermal_conductivity(
                    material_1,
                    temperature_1,
                    pressure_1)?;

                let mut liquid_mcp = mcp_1;
                let mut solid_mcp = mcp_2;

                // we assume liquid is material 1, solid is material 
                // 2 
                // if we are wrong, swop it around
                let mut alpha_liquid: DiffusionCoefficient = 
                thermal_diffusivity(
                    material_1,
                    temperature_1,
                    pressure_1)?;
                
                let mut alpha_solid: DiffusionCoefficient = 
                thermal_diffusivity(
                    material_2,
                    temperature_2,
                    pressure_2)?;

                match material_1 {
                    Solid(_) => {
                        // if material 1 is solid, then we should obtain 
                        // liquid thermal conductivity from material 2
                        liquid_thermal_conductivity = thermal_conductivity(
                            material_2,
                            temperature_2,
                            pressure_2
                        )?;

                        alpha_solid = thermal_diffusivity(
                                material_1,
                                temperature_1,
                                pressure_1)?;

                        alpha_liquid = thermal_diffusivity(
                                material_2,
                                temperature_2,
                                pressure_2)?;

                        liquid_mcp = mcp_2;
                        solid_mcp = mcp_1;
                        ()
                    },
                    Liquid(_) => {
                        // if it's liquid material, match it and obtain 
                        // thermal conductivity 
                        // if the liquid is in material 1, do nothing
                        ()
                    },
                }

                // probably should be taken as cube root but maybe later 
                let liquid_lengthscale: Length = solid_lengthscale * 
                    liquid_mcp / 
                    solid_mcp;


                let solid_timescale: Time = max_mesh_fourier_number * 
                    solid_lengthscale * 
                    solid_lengthscale / 
                    alpha_solid ;


                // Now, I'm going to use the liquid lengthscale 
                // rather than the radius as the lengthscale is a measure 
                // of the thermal inertia 
                //
                // whether it's linear or cube root, I'm not too sure now 
                // didn't bother checking or otherwise. 
                //
                // I'm pretty sure it's cube rooted
                // but I'll correct that later if we get issues

                let nusselt_number: Ratio = 
                h * liquid_lengthscale / liquid_thermal_conductivity;

                let liquid_timescale: Time = max_mesh_fourier_number * 
                    liquid_lengthscale * 
                    liquid_lengthscale / 
                    alpha_liquid / 
                    nusselt_number;

                // now set the control volume using the minimum timescales

                let interaction_minimum_timescale: Time; 
                
                if solid_timescale > liquid_timescale {
                    interaction_minimum_timescale = liquid_timescale;
                } else {
                    interaction_minimum_timescale = solid_timescale;
                }


                // take minimum of the timescales and assign it to 
                // the cv 1 or cv 2 timescales
                //

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }

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
            _cylinder_length)
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

                    let mut error_str: String = "the inner diameter  \n
                        plus shell thicknesses do not equate  \n
                        to outer diameter".to_string();

                    error_str += "supplied outer diameter (m):";
                    error_str += &od.value.to_string();
                    error_str += "expected outer diameter (m):";
                    error_str += &expected_outer_diameter.value.to_string();


                    return Err(error_str);
                }

                // for this, it should be quite straight forward, 
                // find both alphas 
                //
                // by convention, material 1 is inner shell 
                // and material 2 is outer shell 

                let inner_shell_material = inner_material; 
                let inner_shell_temperature = temperature_1;
                let inner_shell_pressure = pressure_1;

                let outer_shell_material = outer_material; 
                let outer_shell_temperature = temperature_2;
                let outer_shell_pressure = pressure_2;


                let alpha_inner: DiffusionCoefficient = 
                thermal_diffusivity(inner_shell_material,
                    inner_shell_temperature,
                    inner_shell_pressure)?;

                let alpha_outer: DiffusionCoefficient = 
                thermal_diffusivity(outer_shell_material,
                    outer_shell_temperature,
                    outer_shell_pressure)?;

                // now let's calculate timescales of both shells 
                // this pertains to conduction only, no convection

                let inner_shell_timescale: Time = max_mesh_fourier_number * 
                inner_thickness * 
                inner_thickness / 
                alpha_inner;

                let outer_shell_timescale: Time = max_mesh_fourier_number * 
                outer_thickness * 
                outer_thickness / 
                alpha_outer;


                // now adjust timestep again

                let interaction_minimum_timescale: Time; 
                
                if outer_shell_timescale > inner_shell_timescale {
                    interaction_minimum_timescale = inner_shell_timescale;
                } else {
                    interaction_minimum_timescale = outer_shell_timescale;
                }

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }


            },
        HeatTransferInteractionType::UserSpecifiedHeatAddition  
            => {
                // user specified heat addition does not make sense 
                // for cv to cv interaction
                return Err("interaction type needs to be \n 
                    thermal conductance".to_string());
            },

        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCustomArea(_) => {

                // user specified heat flux does not make sense 
                // for cv to cv interaction
                return Err("interaction type needs to be \n 
                    thermal conductance".to_string());

            },
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalOuterArea(_,_) => {

                // user specified heat flux does not make sense 
                // for cv to cv interaction
                return Err("interaction type needs to be \n 
                    thermal conductance".to_string());

            },
        HeatTransferInteractionType::
            UserSpecifiedHeatFluxCylindricalInnerArea(_,_) => {

                // user specified heat flux does not make sense 
                // for cv to cv interaction
                return Err("interaction type needs to be \n 
                    thermal conductance".to_string());

            },
        HeatTransferInteractionType::
            DualCartesianThermalConductance(
            (_material_1, thickness_1),
            (_material_2,thickness_2)) => { 

                let thickness_1_length: Length = thickness_1.into();
                let thickness_2_length: Length = thickness_2.into();


                // now, we are not quite sure which material 
                // corresponds to which control volume,
                // though of course, we know both are solids 
                // 
                // and that 
                // Delta t = Fo * Delta x * Delta x  / alpha 
                //
                // so to get the minimum delta t, we take the 
                // shorter of the two lengths, and the larger of 
                // the two alphas 
                //

                let mut larger_alpha: DiffusionCoefficient = 
                alpha_1;

                if alpha_1 < alpha_2 {
                    larger_alpha = alpha_2;
                }


                let mut smaller_lengthscale: Length = thickness_1.into();

                if thickness_2_length < thickness_1_length {
                    smaller_lengthscale = thickness_2_length;
                }

                let interaction_minimum_timescale: Time 
                = max_mesh_fourier_number * 
                smaller_lengthscale *
                smaller_lengthscale / 
                larger_alpha;

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }

            },
        HeatTransferInteractionType::
            DualCartesianThermalConductanceThreeDimension(
            data_dual_cartesian_conduction) 
            => {


                let thickness_1 = 
                data_dual_cartesian_conduction .thickness_1;

                let thickness_2 = 
                data_dual_cartesian_conduction .thickness_2;


                let thickness_1_length: Length = thickness_1.into();
                let thickness_2_length: Length = thickness_2.into();


                // now, we are not quite sure which material 
                // corresponds to which control volume,
                // though of course, we know both are solids 
                // 
                // and that 
                // Delta t = Fo * Delta x * Delta x  / alpha 
                //
                // so to get the minimum delta t, we take the 
                // shorter of the two lengths, and the larger of 
                // the two alphas 
                //

                let mut larger_alpha: DiffusionCoefficient = 
                alpha_1;

                if alpha_1 < alpha_2 {
                    larger_alpha = alpha_2;
                }


                let mut smaller_lengthscale: Length = thickness_1.into();

                if thickness_2_length < thickness_1_length {
                    smaller_lengthscale = thickness_2_length;
                }

                let interaction_minimum_timescale: Time 
                = max_mesh_fourier_number * 
                smaller_lengthscale *
                smaller_lengthscale / 
                larger_alpha;

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }

            },

        HeatTransferInteractionType::
            UserSpecifiedConvectionResistance(
            data_convection_resistance) 
            => {

                // repeat analysis as per user specified thermal 
                // conductance

                let heat_transfer_coeff: HeatTransfer = 
                data_convection_resistance.heat_transfer_coeff;
                let surf_area: Area = 
                data_convection_resistance.surf_area.into();

                let user_specified_conductance: ThermalConductance = 
                heat_transfer_coeff * surf_area;

                let lengthscale_stability_vec_1 = 
                single_cv_1.mesh_stability_lengthscale_vector.clone();
                let lengthscale_stability_vec_2 = 
                single_cv_2.mesh_stability_lengthscale_vector.clone();

                let mut min_lengthscale = Length::new::<meter>(0.0);

                for length in lengthscale_stability_vec_1.iter() {
                    if *length > min_lengthscale {
                        min_lengthscale = *length;
                    }
                }
                
                for length in lengthscale_stability_vec_2.iter() {
                    if *length > min_lengthscale {
                        min_lengthscale = *length;
                    }
                }

                // now that we've gotten both length scales and gotten 
                // the shortest one, we'll need to obtain a timescale 
                // from the conductance 
                //
                // Delta t = Fo * rho * cp * Delta x * resistance
                //
                // we'll convert conductance into total resistance first 
                // 

                let total_resistance = 1.0/user_specified_conductance;

                // then we'll approximate the thermal resistance of 
                // cv 1 using the ratio of heat capacities 
                //
                // This is approximation, not exact, I'll deal with 
                // the approximation as I need to later.

                let approx_thermal_resistance_1= total_resistance *
                heat_capacity_ratio_cv_1;

                let approx_thermal_cond_1: ThermalConductance = 
                1.0/approx_thermal_resistance_1;

                // do note that thermal conductance is in watts per kelvin 
                // not watts per m2/K
                // area is factored in
                //
                // qA = H (T_1 - T_2)
                //
                // if i want HbyA or HbyV, then divide by the control volume
                //
                // H = kA/L for thermal conductivity case
                //  
                // bummer, no surface area!
                // I suppose I'll just use the minimum length scales 
                // then for conservative-ness

                let approx_thermal_conductivity_1: ThermalConductivity = 
                approx_thermal_cond_1 / min_lengthscale;



                let approx_thermal_resistance_2 = total_resistance *
                heat_capacity_ratio_cv_2;

                let approx_thermal_cond_2: ThermalConductance = 
                1.0/approx_thermal_resistance_2;

                let approx_thermal_conductivity_2: ThermalConductivity = 
                approx_thermal_cond_2 / min_lengthscale;

                // let's get timescale 1 and 2 estimate 

                let density_1: MassDensity 
                = density(material_1, temperature_1, pressure_1)?;

                let density_2: MassDensity 
                = density(material_2, temperature_2, pressure_2)?;

                let rho_cp_1: VolumetricHeatCapacity = density_1 * cp_1;
                let rho_cp_2: VolumetricHeatCapacity = density_2 * cp_2;

                let approx_alpha_1: DiffusionCoefficient = 
                approx_thermal_conductivity_1/rho_cp_1;

                let approx_alpha_2: DiffusionCoefficient = 
                approx_thermal_conductivity_2/rho_cp_2;


                // Fo  = alpha * Delta t / Delta x / Delta x 
                //
                // Delta t = Fo * Delta x * Delta x / alpha 
                //

                let timescale_1: Time = max_mesh_fourier_number 
                * min_lengthscale 
                * min_lengthscale 
                / approx_alpha_1;

                let timescale_2: Time = max_mesh_fourier_number 
                * min_lengthscale 
                * min_lengthscale 
                / approx_alpha_2;

                let interaction_minimum_timescale: Time; 
                
                if timescale_1 > timescale_2 {
                    interaction_minimum_timescale = timescale_2;
                } else {
                    interaction_minimum_timescale = timescale_1;
                }


                // take minimum of the timescales and assign it to 
                // the cv 1 or cv 2 timescales
                //

                if cv_1_timestep > interaction_minimum_timescale {
                    cv_1_timestep = interaction_minimum_timescale;
                }

                if cv_2_timestep > interaction_minimum_timescale {
                    cv_2_timestep = interaction_minimum_timescale;
                }


            },
    };

    // push the corrected minimum timesteps to cv 1 and cv 2
    single_cv_1.max_timestep_vector.push(cv_1_timestep);
    single_cv_2.max_timestep_vector.push(cv_2_timestep);

    // load the minimum timestep and return
    let mut minimum_timestep: Time = cv_1_timestep;

    if cv_1_timestep > cv_2_timestep {
        minimum_timestep = cv_2_timestep;
    }

    return Ok(minimum_timestep);
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
                // if the control volume is fluid, we will need 
                // to introduce another time scale

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
                // if the control volume is fluid, we will need 
                // to introduce another time scale

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
