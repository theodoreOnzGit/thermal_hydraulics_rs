use uom::si::length::meter;
use uom::{si::ratio::ratio, ConstZero};
use uom::si::f64::*;

use super::{enums::NusseltCorrelation, input_structs::GnielinskiData};

/// from Du's paper
///
/// Du, B. C., He, Y. L., Qiu, Y., Liang, Q., & Zhou, Y. P. (2018). 
/// Investigation on heat transfer characteristics of molten salt in 
/// a shell-and-tube heat exchanger. International Communications 
/// in Heat and Mass Transfer, 96, 61-68.
///
/// we have a generic Gnielinski type correlation, 
/// empirically fitted to experimental data. This is in the form:
///
/// Nu = C (Re^m - 280.0) Pr_f^0.4 ( 1.0 + (D_e/l)^(2/3) ) ( Pr_f / Pr_w )^0.25
///
/// For Du's Heat exchanger, 
/// C = 0.04318,
/// m = 0.7797
///
/// Re, Nu_shell
/// 3510.033,42.582
/// 3571.349,43.32
/// 3691.75,43.852
/// 3751.951,44.672
/// 3794.314,44.795
/// 3847.826,45.574
/// 3959.309,47.09
/// 4019.509,47.459
/// 4267.001,53.238
/// 4356.187,54.836
/// 4550.167,58.238
/// 4630.435,59.303
/// 4730.769,60.451
/// 4942.586,62.582
/// 5230.212,63.77
/// 5388.517,64.344
/// 5481.048,65.861
///
/// From the paper 
/// for the salt, temperatures range from 204-236 C 
/// Pr is from 19.82 to 24.03, 
///
/// Pr = 22 seems reasonable for bulk fluid (Pr_f)
///
/// and the correction factor Pr_f/Pr_w is from 
/// 0.4273 to 0.5646, decent estimate is 0.5
///
/// These values allow us to calculate the salt Nusselt numbers 
/// to reproduce the Re and Nu_shell data.
///
/// I'm going to try the values at Re = 3510 (which is in the transitional 
/// regime), Re = 4019, which is in turbulent regime, 
/// and Re = 5481, which is also in the turbulent regime
///
///
#[test] 
pub fn du_correlation_empirical_test(){

    let c = Ratio::new::<ratio>(0.04318);
    let m = 0.7797_f64;

    // now some parameters to determine things,
    // from Du's paper
    
    let heat_exchg_length = Length::new::<meter>(1.95);
    let tube_od = Length::new::<meter>(0.014);
    let shell_id = Length::new::<meter>(0.1);

    let number_of_tubes = 19;

    // from Du's paper, eqn 14
    // D_e = (D_i^2 - N_t d_o^2)/(D_i + N_t d_o)
    let effective_diameter = (
        shell_id * shell_id - 
        number_of_tubes as f64 
        * tube_od * tube_od) / 
        (shell_id + number_of_tubes as f64 * tube_od);
    let length_to_diameter = heat_exchg_length/effective_diameter;

    let gnielinski_params: GnielinskiData = 
        GnielinskiData { 
            reynolds: Ratio::ZERO, 
            prandtl_bulk: Ratio::ZERO, 
            prandtl_wall: Ratio::ZERO, 
            darcy_friction_factor: Ratio::ZERO, 
            length_to_diameter,
        };



    let du_nusselt_correlation: NusseltCorrelation 
        = NusseltCorrelation::CustomGnielinskiGeneric(
            gnielinski_params, c, m);

    let bulk_prandtl_number = Ratio::new::<ratio>(22.0);
    // Pr_f/Pr_w  is about 0.5, hence 
    // Pr_w/Pr_f = 2.0 approx 
    //
    // Pr_w = Pr_w/Pr_f * Pr_f
    let wall_prandtl_number = bulk_prandtl_number / 0.5;


    // test for Re about 5481
    {
        let reynolds = Ratio::new::<ratio>(5481_f64);
        let nusselt = du_nusselt_correlation
            .estimate_based_on_prandtl_reynolds_and_wall_correction(
                bulk_prandtl_number, 
                wall_prandtl_number, 
                reynolds).unwrap();
    
        approx::assert_relative_eq!(
            nusselt.get::<ratio>(),
            0.0,
            max_relative=0.00
            );
    
    }



}
