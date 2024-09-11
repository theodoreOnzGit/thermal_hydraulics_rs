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
/// on the shell side.
/// There is parasitic heat loss through the heater when there should be 
/// none. (See Zou's publication on nuclear engineering and design in 2021)
/// this serves as a baseline as to what kind of heat losses to expect
pub mod coupled_dracs_loop_ver_1;

/// for the coupled dracs loop, we need to calibrate heat loss 
/// through the primary loop 
/// hot leg (from heater outlet to dhx shell inlet)
/// and cold leg 
/// (from dhx shell outlet to heater inlet)
/// 
/// The data in csv format from Zweibaum's unpublished work in Th Lab 
/// Archives is:
///
/// dataset number,pri loop mass flowrate (kg/s),Heater outlet (DegC),DHX shell top (DegC),DHX shell bottom (DegC),Heater inlet (DegC),
/// C-1,0.02003,75.22747,71.47752,53.60943,50.45784,
/// C-2,0.02367,82.41863,78.36713,57.13467,53.79036,
/// C-3,0.02635,87.78188,84.37342,59.82845,56.71891,
/// C-4,0.02949,94.71628,90.97595,63.9812,60.83029,
/// C-5,0.0319,100.37023,96.20228,67.05336,64.07406,
/// C-6,0.03412,105.25073,101.3375,69.85085,67.1654,
/// C-7,0.03562,110.34289,106.43149,73.21226,70.6215,
/// C-8,0.03593,115.52364,111.37615,76.13202,73.63344,
/// C-9,0.03547,119.96879,116.05003,79.02407,76.54479,
///
/// The basic procedure is to fix the mass flowrate through the branches
/// and then use a PID controller to set a heat transfer coefficient  
/// so that the inlet dhx shell temperature is equal to the set point,
/// That of the experimental data
///
/// repeat the same for the cold leg.
///
pub mod primary_loop_parasitic_heat_loss_calibration;
