
#[macro_use]
extern crate approx;
use uom::si::f64::*;

#[test]
pub fn meant_to_fail(){

    // now for rust , we don't have assert equal
    // showing expected and test values
    // we just see if left == right
    // not like C#,
    // where left is expected value,
    // right is asserted value
    //
    assert_eq!(2.0,2.0);
    unimplemented!();
}

/// This test prototypes the CIET 
#[test]
pub fn ciet_crude_heater_v_1_0 (){

    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::mass_rate::kilogram_per_second;
    // for each timestep, the ciet heater must have 
    // the inlet conditions specified
    let fluid_temp_inlet = ThermodynamicTemperature::new::<
        degree_celsius>(79.0);
    let mass_flowrate = MassRate::new::<kilogram_per_second>(0.18);

    // we also need to let it know that therminol VP1 is the correct 
    // fluid 
    //
    // This can be put into some enum kind of thing, because 
    // writing those correlations is rather repetitive
    //
    // suppose that is done, we need to use this to specify the 
    // inlet enthalpy 
    //
    // at the bare minimum I need a function to convert the temperature 
    // to an enthalpy
    // 
    // this is found in the fluid_mechanics_lib for now 

    use thermal_hydraulics_rs::fluid_mechanics_lib::prelude::*;
    use thermal_hydraulics_rs::fluid_mechanics_lib::therminol_component;

    // we'll get the inlet specific enthalpy 

    let inlet_specific_enthalpy = therminol_component::dowtherm_a_properties:: 
        getDowthermAEnthalpy(fluid_temp_inlet);

    // now we calculate inlet enthalpy
    use thermal_hydraulics_rs::heat_transfer_lib:: 
    control_volume_calculations::common_functions;

    let inlet_enthalpy = common_functions::calculate_enthalpy_flow(
        mass_flowrate, inlet_specific_enthalpy);

    // print some output
    println!("{:?}",inlet_enthalpy);

    // next step, need to calculate heat flow between fluid and 
    // environment, so we need a fluid temperature at 
    // present timestep 
    //
    // Also need a heater_shell temperature at present timestep
    //
    let fluid_temperature_present_timestep = ThermodynamicTemperature:: 
        new::<degree_celsius>(88.0);

    let heater_shell_temperature_present_timestep = 
    ThermodynamicTemperature::new::<degree_celsius>(138.0);

    // we can assume the heater shell is a lumped capacitance model 
    // or something
    // otherwise there should be some temperature gradient between 
    // heater center, heater inner surface and heater outer surface
    //

}
