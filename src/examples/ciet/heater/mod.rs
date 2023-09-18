//! The heater module here represents components from 
//!
//! BT-11 to BT-12
//! This is because heater inlet and outlet temperatures are measured 
//! using BT-11 and BT-12
//! 
//! However, BT-11 is located before the heater bottom head 
//! and BT-12 is located after the static mixer MX-10 
//!
//! Hence, we have four sections,
//!
//! the heater top head, heater bottom head, heated section,
//! mixer MX-10 and the static mixer pipe attached to MX-10 modelled 
//! in the Relap and SAM model
//!
//! So there is not only some residence time, but also other mechanisms 
//! for parasitic heat loss 
//!
//! Dr Dane De Wet's transform model does callibrate for this using 
//! a heat transfer coefficient of 20 W/(m^2 K) instead of 6 W/(m^2 K)
//!
//!
//! I intend to connect structural supports to the heater top and bottom 
//! head and callibrate the length of those structural supports 
//! as part of model callibration such that the heater inlet is 
//! 80C approximately, and the heater outlet is 102.45 C 
//!
//! at nominal heater power of 8 kW
//!
//! For this, I also want to ensure that the code runs fast enough,
//! at least faster than real time, so it is suitable for digital 
//! twin applications
//!
//!
//! 
//!
/// represents heater version 1, 
///
/// it allows fluid to flow through it in an annular tube
pub struct HeaterVersion1;

/// represents heater version 2, 
///
/// it has twisted tape
pub struct HeaterVersion2;

pub mod heater_version_2_bare;
pub use heater_version_2_bare::*;

pub mod heater_top_and_bottom_head_bare;
pub use heater_top_and_bottom_head_bare::*;

pub mod static_mixer_mx_10;
pub use static_mixer_mx_10::*;


pub fn example_heater(){

    // construct

    let _heater_v1 = HeaterVersion1{};
    let _heater_v2 = HeaterVersion2{};
    //let _heater_v2_bare = HeaterVersion2Bare{};

}
#[test]
pub fn example_heater_test(){

}
