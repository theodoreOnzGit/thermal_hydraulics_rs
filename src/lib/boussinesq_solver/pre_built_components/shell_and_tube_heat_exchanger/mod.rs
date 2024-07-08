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
