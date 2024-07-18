// I'm taking data from:
// Du, Bao-Cun, et al. "Investigation on heat transfer 
// characteristics of molten salt in a shell-and-tube heat 
// exchanger." International Communications in Heat and 
// Mass Transfer 96 (2018): 61-68.
//
// to validate the models here. Hopefully it works!
//
// Basically, I'm not given data straight up, nor am I given the conditions 
// used to produce the data. However, I am given some data points of 
// the external tube nusselt number, which correspond to experimentally 
// determined nusselt numbers for the salt.
//
// I'm also given a Nusselt correlation for the shell side 
// (outside the tube) which contains the salt 
//
// Nu = 0.04318 * (Re^0.7797 - 280) Pr ^(0.4) ( 1 + (D_e/l)^(2/3) ) * (Pr_f/Pr_w)^0.25
//
// If i use this correlation for the outside of the tube, I should technically 
// be able to reproduce the results.
//
// The Re and Nu_shell numbers (Fig 9, obtained using graph reader) are:
//
// Re, Nu_shell
// 3510.033,42.582
// 3571.349,43.32
// 3691.75,43.852
// 3751.951,44.672
// 3794.314,44.795
// 3847.826,45.574
// 3959.309,47.09
// 4019.509,47.459
// 4267.001,53.238
// 4356.187,54.836
// 4550.167,58.238
// 4630.435,59.303
// 4730.769,60.451
// 4942.586,62.582
// 5230.212,63.77
// 5388.517,64.344
// 5481.048,65.861
//
// So basically, we need to reproduce this data. Use the correlations 
// and obtain the Nu_shell using some of the experiments. This is a 
// verification study.
//
//
// I'll still need to program in the Nusselt correlations, as well 
// as the as the thermophysical properties 
//
// Also, what do I use as the input parameters?
//
// I'm only given a flow rate of molten salt:
// 12.63 - 16.23 m3/h
//
// Flow rate of oil (almost constant):
// 15.48 - 15.79 m3/h
// 
// inlet temperature of molten salt (HITEC):
// 214.93 - 236.91 (deg C)
//
// outlet temperature of molten salt (HITEC):
// 204.06 - 225.35 (deg C)
//
// Inlet Temperature of oil:
// 74.49 - 90.41 (deg C)
//
// Outlet Temperature of oil:
// 88.24 - 110.74 (deg C) 
//
// And the Reynold's number on the shell side is meant to be:
// Re_s = rho_s u_s D_e / manual_tests
// 
// The superifical velocity is:
//
// u_s = q_(m,s) / (rho_s A_c)
//
// The cross sectional area is the total shell inner area
// minus the area of the area of the tubes:
//
// A_c = pi/4 * D_i^2 - N_t * pi/4 * d_o^2
//
// N_t = number of tubes
//
// The effective hydraulic diameter for the tube is:
//
// D_e = 4 * A_c / P_w
//
// where P_w  is the wetted perimeter
//
// The wetted (heated) perimeter is: 
// Pi * D_i + N_t * Pi * d_o
//
// Substitution results in a final expression:
//
// D_e = (D_i^2 - N_t d_o^2) / (D_i + N_t d_o)
//
// The Re for the oil is from 3514 - 5482, and Pr is from 19.82 to 24.03
//
// For programming, we need two sets of fluid entities within this 
// heat exchanger.
//
// For the reference simulation work, the tube is assumed to be well 
// insulated. (adiabatic BC)
//
// There are a few steps to complete to try and replicate the experiment
//


use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;

use super::heat_transfer_entities::cv_types::CVType;
use super::heat_transfer_entities::HeatTransferEntity;
use uom::si::f64::*;

/// Single pass, no baffle, parrallel flow shell and tube 
/// heat exchanger 
///
/// Nusselt correlations can be customised to empirically fit other 
/// correlations
///
/// The axial sides are adiabatic unless otherwise stated
#[derive(Clone,Debug,PartialEq)]
pub struct SimpleShellAndTubeHeatExchanger {

    
    inner_nodes: usize,

    /// this HeatTransferEntity represents the pipe shell which 
    /// contains the tube side fluid
    /// it is exposed to the shell side fluid
    pub pipe_shell: HeatTransferEntity,


    /// this HeatTransferEntity represents the pipe fluid
    /// which is coupled to the pipe shell via a Nusselt Number based
    /// thermal resistance (usually Gnielinski correlation)
    pub tube_side_parallel_fluid_array: HeatTransferEntity,

    /// this HeatTransferEntity represents the pipe fluid
    /// which is coupled to the pipe shell via a Nusselt Number based
    /// thermal resistance, this must be specified by the user
    pub shell_side_fluid_array: HeatTransferEntity,

    /// this HeatTransferEntity represents the pipe shell which 
    /// contains the shell side fluid and tube bundle
    /// it is exposed to the insulation, or ambient temperature 
    /// depending on whether the insulation is toggled on or off
    pub outer_shell: HeatTransferEntity,

    /// ambient temperature that the shell and tube heat 
    /// exchanger is exposed to
    pub ambient_temperature: ThermodynamicTemperature,

    /// heat transfer coefficient to ambient
    /// This provides thermal resistance between the surface of 
    /// the shell and tube heat exchanger 
    /// This could be the outer shell or insulation, depending on whether 
    /// insulation is toggled on or off
    /// 
    ///
    pub heat_transfer_to_ambient: HeatTransfer,

    /// insulation array covering the 
    /// outer_shell array if insulation is toggled on
    pub insulation_array: HeatTransferEntity,

    /// this option allows the user to toggle on or off insulation 
    pub heat_exchanger_has_insulation: bool,

    /// representative 
    /// tube outer diameter on a per tube bases
    pub tube_side_od: Length,

    /// representative 
    /// tube inner diameter one a per tube basis
    pub tube_side_id: Length,

    /// representative tube flow area on a per tube basis
    pub tube_side_flow_area: Area,

    /// loss correlation on a per tube basis
    pub tube_side_custom_component_loss_correlation: DimensionlessDarcyLossCorrelations,

    /// number of tubes in parallel 
    /// each pipe fluid array represents one tube only
    pub number_of_tubes: u32,

    /// assuming the outer shell is circular, provide the internal diameter 
    pub shell_side_id: Length,

    /// assuming the outer shell is circular, provide the outer diameter 
    pub shell_side_od: Length,

    /// allows for a custom flow area for the shell side
    pub shell_side_flow_area: Area,


    /// specifies an thickness for the insulation covering 
    /// the shell side
    pub insulation_thickness: Length,


}


/// stuff such as conductances are calculated here
pub mod preprocessing;

/// implementations for the FluidComponent trait
/// are done here
pub mod fluid_component;


/// stuff for calculation is done here, ie, advancing timestep
pub mod calculation;

/// postprocessing stuff, ie, get the temperature vectors 
/// of both arrays of control volumes 
pub mod postprocessing;

/// type conversion, such as into fluid component and such
pub mod type_conversion;


/// verification tests for parallel tubing
pub mod tests;
