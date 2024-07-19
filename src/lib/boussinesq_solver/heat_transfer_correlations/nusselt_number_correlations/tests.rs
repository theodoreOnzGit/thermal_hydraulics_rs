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
/// 0.4273 to 0.5646
///
/// These values allow us to calculate the salt Nusselt numbers 
/// to reproduce the Re and Nu_shell data.
///
///
#[test] 
pub fn du_correlation_empirical_test(){

}
