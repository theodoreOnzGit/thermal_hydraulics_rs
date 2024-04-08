use std::f64::consts::PI;

use uom::si::angle::degree;
use uom::si::area::square_meter;
// This library was developed for use in my PhD thesis under supervision 
// of Professor Per F. Peterson. It is part of a thermal hydraulics
// library in Rust that is released under the GNU General Public License
// v 3.0. This is partly due to the fact that some of the libraries 
// inherit from GeN-Foam and OpenFOAM, both licensed under GNU General
// Public License v3.0.
//
// As such, the entire library is released under GNU GPL v3.0. It is a strong 
// copyleft license which means you cannot use it in proprietary software.
//
//
// License
//    This file is part of thermal_hydraulics_rs, a partial library of the
//    thermal hydraulics library written in rust meant to help with the
//    fluid mechanics aspects of the calculations
//     
//    Copyright (C) 2022-2023  Theodore Kay Chen Ong, Singapore Nuclear
//    Research and Safety Initiative, Per F. Peterson, University of 
//    California, Berkeley Thermal Hydraulics Laboratory
//
//    thermal_hydraulics_rs is free software; you can redistribute it and/or modify it
//    under the terms of the GNU General Public License as published by the
//    Free Software Foundation; either version 2 of the License, or (at your
//    option) any later version.
//
//    thermal_hydraulics_rs is distributed in the hope that it will be useful, but WITHOUT
//    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//    for more details.
//
//    This library is part of a thermal hydraulics library in rust
//    and contains some code copied from GeN-Foam, and OpenFOAM derivative.
//    This offering is not approved or endorsed by the OpenFOAM Foundation nor
//    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
//    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.
//    Nor is it endorsed by the authors and owners of GeN-Foam.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// © All rights reserved. Theodore Kay Chen Ong,
// Singapore Nuclear Research and Safety Initiative,
// Per F. Peterson,
// University of California, Berkeley Thermal Hydraulics Laboratory
//
// Main author of the code: Theodore Kay Chen Ong, supervised by
// Professor Per F. Peterson
use uom::si::f64::*;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::length::{meter, millimeter};
use uom::si::ratio::ratio;
use uom::si::thermodynamic_temperature::degree_celsius;
use uom::si::pressure::atmosphere;

use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;

use super::insulated_pipes_and_fluid_components::InsulatedFluidComponent;
use super::non_insulated_fluid_components::NonInsulatedFluidComponent;

/// creates a new pipe6a for CIET using the RELAP5-3D and SAM parameters 
/// Pipe6a in Compact Integral Effects Test (CIET)
/// CTAH branch 
///
/// It is a static mixer pipe
/// otherwise known as the static mixer pipe 6a
///
/// Zou, Ling, Rui Hu, and Anne Charpentier. SAM code 
/// validation using the compact integral effects test (CIET) 
/// experimental data. No. ANL/NSE-19/11. Argonne National Lab.(ANL), 
///
///
/// Zweibaum, Nicolas. Experimental validation of passive safety 
/// system models: Application to design and optimization of 
/// fluoride-salt-cooled, high-temperature reactors. University of 
/// California, Berkeley, 2015.
/// Argonne, IL (United States), 2019.
pub fn new_pipe_6a() -> InsulatedFluidComponent {
    let initial_temperature = ThermodynamicTemperature::new::<degree_celsius>(21.7);
    let ambient_temperature = ThermodynamicTemperature::new::<degree_celsius>(20.0);
    let fluid_pressure = Pressure::new::<atmosphere>(1.0);
    let solid_pressure = Pressure::new::<atmosphere>(1.0);
    let hydraulic_diameter = Length::new::<meter>(2.79e-2);
    let pipe_length = Length::new::<meter>(0.1526);
    let flow_area = hydraulic_diameter * hydraulic_diameter * PI/4.0;
    let incline_angle = Angle::new::<degree>(51.526384);
    let form_loss = Ratio::new::<ratio>(5.05);
    //estimated component wall roughness (doesn't matter here,
    //but i need to fill in)
    let surface_roughness = Length::new::<millimeter>(0.015);
    let shell_id = hydraulic_diameter;
    let pipe_thickness = Length::new::<meter>(0.0027686);
    let shell_od = shell_id + pipe_thickness;
    let insulation_thickness = Length::new::<meter>(0.0508);
    let pipe_shell_material = SolidMaterial::SteelSS304L;
    let insulation_material = SolidMaterial::Fiberglass;
    let pipe_fluid = LiquidMaterial::TherminolVP1;
    let htc_to_ambient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // from SAM nodalisation, we have 2 nodes only, 
    // now because there are two outer nodes, the 
    // number of inner nodes is zero
    let user_specified_inner_nodes = 0; 

    let insulated_component = InsulatedFluidComponent::new_insulated_pipe(
        initial_temperature, 
        ambient_temperature, 
        fluid_pressure, 
        solid_pressure, 
        flow_area, 
        incline_angle, 
        form_loss, 
        shell_id, 
        shell_od, 
        insulation_thickness, 
        pipe_length, 
        hydraulic_diameter, 
        pipe_shell_material, 
        insulation_material, 
        pipe_fluid, 
        htc_to_ambient, 
        user_specified_inner_nodes, 
        surface_roughness);

    insulated_component
}

/// creates a new pipe6a for CIET using the RELAP5-3D and SAM parameters 
/// Component 6 in Compact Integral Effects Test (CIET)
/// CTAH branch  (also known as static mixer 41)
///
/// static mixer 41
/// label component 6 
/// in Compact Integral Effects Test (CIET)
/// CTAH branch 
/// static mixer 41 (MX-41) on CIET diagram
/// in the pump and CTAH branch
/// just before CTAH (AKA IHX)
/// from top to bottom
///
/// label 6 on diagram
///
/// It is a static mixer pipe
/// otherwise known as the static mixer pipe 6a
///
/// Zou, Ling, Rui Hu, and Anne Charpentier. SAM code 
/// validation using the compact integral effects test (CIET) 
/// experimental data. No. ANL/NSE-19/11. Argonne National Lab.(ANL), 
///
///
/// Zweibaum, Nicolas. Experimental validation of passive safety 
/// system models: Application to design and optimization of 
/// fluoride-salt-cooled, high-temperature reactors. University of 
/// California, Berkeley, 2015.
/// Argonne, IL (United States), 2019.
pub fn new_static_mixer_41() -> InsulatedFluidComponent {
    let initial_temperature = ThermodynamicTemperature::new::<degree_celsius>(21.7);
    let ambient_temperature = ThermodynamicTemperature::new::<degree_celsius>(20.0);
    let fluid_pressure = Pressure::new::<atmosphere>(1.0);
    let solid_pressure = Pressure::new::<atmosphere>(1.0);
    let hydraulic_diameter = Length::new::<meter>(2.79e-2);
    let component_length = Length::new::<meter>(0.33);
    let flow_area = Area::new::<square_meter>(6.11e-4);
    let incline_angle = Angle::new::<degree>(51.526384);
    let form_loss = Ratio::new::<ratio>(21.0);
    let reynolds_power = -1_f64;
    let reynolds_coefficient = Ratio::new::<ratio>(4000.0);
    //estimated component wall roughness (doesn't matter here,
    //but i need to fill in)
    let shell_id = hydraulic_diameter;
    let pipe_thickness = Length::new::<meter>(0.0027686);
    let shell_od = shell_id + pipe_thickness;
    let insulation_thickness = Length::new::<meter>(0.0508);
    let pipe_shell_material = SolidMaterial::SteelSS304L;
    let insulation_material = SolidMaterial::Fiberglass;
    let pipe_fluid = LiquidMaterial::TherminolVP1;
    let htc_to_ambient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // from SAM nodalisation, we have 2 nodes only, 
    // now because there are two outer nodes, the 
    // number of inner nodes is zero
    let user_specified_inner_nodes = 0; 

    let insulated_component = InsulatedFluidComponent::new_custom_component(
        initial_temperature, 
        ambient_temperature, 
        fluid_pressure, 
        solid_pressure, 
        flow_area, 
        incline_angle, 
        form_loss, 
        reynolds_coefficient, 
        reynolds_power, 
        shell_id, 
        shell_od, 
        insulation_thickness, 
        component_length, 
        hydraulic_diameter, 
        pipe_shell_material, 
        insulation_material, 
        pipe_fluid, 
        htc_to_ambient, 
        user_specified_inner_nodes);

    insulated_component
}

/// creates a new ctah vertical for CIET using the RELAP5-3D and SAM parameters 
/// in Compact Integral Effects Test (CIET)
///
/// this is inactive, so it behaves more like a pipe rather than a 
/// heat exchanger
///
/// Vertical part of Coiled Tube Air Heater (CTAH)
/// label component 7a
/// in Compact Integral Effects Test (CIET)
/// CTAH branch 
///
/// It is NOT insulated by the way
///
/// It is a static mixer pipe
/// otherwise known as the static mixer pipe 6a
///
/// Zou, Ling, Rui Hu, and Anne Charpentier. SAM code 
/// validation using the compact integral effects test (CIET) 
/// experimental data. No. ANL/NSE-19/11. Argonne National Lab.(ANL), 
///
///
/// Zweibaum, Nicolas. Experimental validation of passive safety 
/// system models: Application to design and optimization of 
/// fluoride-salt-cooled, high-temperature reactors. University of 
/// California, Berkeley, 2015.
/// Argonne, IL (United States), 2019.
///
pub fn new_inactive_ctah_vertical() -> NonInsulatedFluidComponent {
    let initial_temperature = ThermodynamicTemperature::new::<degree_celsius>(21.7);
    let ambient_temperature = ThermodynamicTemperature::new::<degree_celsius>(20.0);
    let fluid_pressure = Pressure::new::<atmosphere>(1.0);
    let solid_pressure = Pressure::new::<atmosphere>(1.0);
    let hydraulic_diameter = Length::new::<meter>(1.19e-2);
    let pipe_length = Length::new::<meter>(0.3302);
    let flow_area = Area::new::<square_meter>(1.33E-03);
    let incline_angle = Angle::new::<degree>(-90.0);
    let form_loss = Ratio::new::<ratio>(3.9);
    //estimated component wall roughness (doesn't matter here,
    //but i need to fill in)
    let surface_roughness = Length::new::<millimeter>(0.015);
    let id = hydraulic_diameter;
    let pipe_thickness = Length::new::<meter>(0.000406);
    let od = id + pipe_thickness;
    let pipe_shell_material = SolidMaterial::SteelSS304L;
    let pipe_fluid = LiquidMaterial::TherminolVP1;
    let htc_to_ambient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // from SAM nodalisation, we have 3 nodes only, 
    // now because there are two outer nodes, the 
    // number of inner nodes is zero
    let user_specified_inner_nodes = 3-2; 

    let non_insulated_component = NonInsulatedFluidComponent::new_bare_pipe(
        initial_temperature, 
        ambient_temperature, 
        fluid_pressure, 
        solid_pressure, 
        flow_area, 
        incline_angle, 
        form_loss, 
        id, 
        od, 
        pipe_length, 
        hydraulic_diameter, 
        surface_roughness, 
        pipe_shell_material, 
        pipe_fluid, 
        htc_to_ambient, 
        user_specified_inner_nodes);

    non_insulated_component
}



/// creates a new ctah vertical for CIET using the RELAP5-3D and SAM parameters 
/// in Compact Integral Effects Test (CIET)
///
/// this is inactive, so it behaves more like a pipe rather than a 
/// heat exchanger
///
/// Horizontal part of Coiled Tube Air Heater (CTAH)
/// label component 7b
/// in Compact Integral Effects Test (CIET)
/// CTAH branch 
/// coiled tube air heater
/// has fldk = 400 + 52,000/Re
///
/// label is 7b
/// empirical data in page 48 on pdf viewer in Dr
/// Zweibaum thesis shows reverse flow has same
/// pressure drop characteristics as forward flow
///
/// It is NOT insulated by the way
///
/// It is a static mixer pipe
/// otherwise known as the static mixer pipe 6a
///
/// Zou, Ling, Rui Hu, and Anne Charpentier. SAM code 
/// validation using the compact integral effects test (CIET) 
/// experimental data. No. ANL/NSE-19/11. Argonne National Lab.(ANL), 
///
///
/// Zweibaum, Nicolas. Experimental validation of passive safety 
/// system models: Application to design and optimization of 
/// fluoride-salt-cooled, high-temperature reactors. University of 
/// California, Berkeley, 2015.
/// Argonne, IL (United States), 2019.
///
pub fn new_inactive_ctah_horizontal() -> NonInsulatedFluidComponent {
    let initial_temperature = ThermodynamicTemperature::new::<degree_celsius>(21.7);
    let ambient_temperature = ThermodynamicTemperature::new::<degree_celsius>(20.0);
    let fluid_pressure = Pressure::new::<atmosphere>(1.0);
    let solid_pressure = Pressure::new::<atmosphere>(1.0);
    let hydraulic_diameter = Length::new::<meter>(1.19e-2);
    let component_length = Length::new::<meter>(1.2342);
    let flow_area = Area::new::<square_meter>(1.33E-03);
    let incline_angle = Angle::new::<degree>(0.0);
    let form_loss = Ratio::new::<ratio>(400.0);
    let reynolds_power = -1_f64;
    let reynolds_coefficient = Ratio::new::<ratio>(52000_f64);
    //estimated component wall roughness (doesn't matter here,
    //but i need to fill in)
    let shell_id = hydraulic_diameter;
    let pipe_thickness = Length::new::<meter>(0.000406);
    let shell_od = shell_id + pipe_thickness;
    let pipe_shell_material = SolidMaterial::SteelSS304L;
    let pipe_fluid = LiquidMaterial::TherminolVP1;
    let htc_to_ambient = HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);
    // from SAM nodalisation, we have 11 nodes only, 
    // now because there are two outer nodes, the 
    // number of inner nodes is zero
    let user_specified_inner_nodes = 11-2; 

    let non_insulated_component = NonInsulatedFluidComponent::
        new_custom_component(
            initial_temperature, 
            ambient_temperature, 
            fluid_pressure, 
            solid_pressure, 
            flow_area, 
            incline_angle, 
            form_loss, 
            reynolds_coefficient, 
            reynolds_power, 
            shell_id, 
            shell_od, 
            component_length, 
            hydraulic_diameter, 
            pipe_shell_material, 
            pipe_fluid, 
            htc_to_ambient, 
            user_specified_inner_nodes);

    non_insulated_component
}
