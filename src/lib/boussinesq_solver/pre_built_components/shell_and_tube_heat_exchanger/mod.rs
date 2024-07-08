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
// The Re and Nu_shell numbers (Fig 9) are:
//
// (TBD)
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
// In other words, we need two sets of fluid entities within this 
// heat exchanger
//
// 
//
//
//
// There are a few steps to complete,
//
