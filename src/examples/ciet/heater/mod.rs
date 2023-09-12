/// represents heater version 1, 
///
/// it allows fluid to flow through it in an annular tube
pub struct HeaterVersion1;

/// represents heater version 2, 
///
/// it has twisted tape
pub struct HeaterVersion2;

/// represents heater version 2 without insulation 
/// This is because during 2018-ish, the heater insulation 
/// got burnt off and a lot of frequency response tests were done 
/// with insulation removed
pub struct HeaterVersion2Bare;

//use thermal_hydraulics_rs::prelude::alpha_nightly::*;

pub fn example_heater(){

    // construct

    let _heater_v1 = HeaterVersion1{};
    let _heater_v2 = HeaterVersion2{};
    let _heater_v2_bare = HeaterVersion2Bare{};

}
