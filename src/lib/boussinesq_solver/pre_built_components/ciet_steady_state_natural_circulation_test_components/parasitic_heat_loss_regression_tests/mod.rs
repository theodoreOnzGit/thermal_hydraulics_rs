/// for gnielinski type correlations, there is a wall correction 
/// factor which is something like (Pr_f/Pr_w)^0.11 
///
/// For cooling situations, this lowers the heat transfer 
/// coefficient. To test if this correction factor is working,
/// I compare the parasitic heat loss of the DRACS loop 
/// without correction factors as a base case,
/// then switch it on to see if there is less parasitic 
/// heat loss
///
/// 
pub mod wall_correction_isolated_dracs_loop_regression;

/// version 1 of coupled DRACS loop 
/// for version 1 of coupled DRACS loop
///
/// no calibration is done. heat exchanger correlations are Gnielinksi type 
/// this serves as a baseline as to what kind of heat losses to expect
pub mod coupled_dracs_loop_ver_1;
