/// property correlations for dowtherm_a,
/// also known as therminol vp1
pub mod dowtherm_a;


/// property correlations for hitec (a nitrate salt)
/// Du, Bao-Cun, et al. "Investigation on 
/// heat transfer characteristics of molten salt in a 
/// shell-and-tube heat exchanger." International Communications 
/// in Heat and Mass Transfer 96 (2018): 61-68.
pub mod hitec_nitrate_salt;


/// Thermophysical of the equimolar mixture 
/// NaNO3-KNO3 from 300 to 600°C
///
/// Nissen, Donald A. "Thermophysical properties 
/// of the equimolar mixture sodium nitrate-
/// potassium nitrate from 300 to 600. degree. C
/// ." Journal of Chemical and Engineering 
/// Data 27.3 (1982): 269-273.
///
/// used in:
///
/// Srivastava, A. K., Kudariyawar, J. Y., 
/// Borgohain, A., Jana, S. S., Maheshwari, N. K., & 
/// Vijayan, P. K. (2016). Experimental and 
/// theoretical studies on the natural circulation 
/// behavior of molten salt loop. 
/// Applied Thermal Engineering, 98, 513-521.
///
/// I used the correlation from Srivastava et al. 
pub mod sodium_potassium_equimolar_nissen_nitrate_salt;

/// properties for a custom liquid material 
/// not covered in the database
/// 
/// You'll need to define your own functions for this to work
pub mod custom_liquid_material;
