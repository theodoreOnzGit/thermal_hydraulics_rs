//! The heater here
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


pub fn example_heater(){

    // construct

    let _heater_v1 = HeaterVersion1{};
    let _heater_v2 = HeaterVersion2{};
    //let _heater_v2_bare = HeaterVersion2Bare{};

}
#[test]
pub fn example_heater_test(){

}
