
#[macro_use]
extern crate approx;
use std::f32::consts::PI;

use uom::si::{f64::*, heat_transfer::watt_per_square_meter_kelvin};

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
    // We could simulate this using some kind of node system, 
    // and we have then two finite volumes each with a lumped capacitance
    // system

    let heater_outer_shell_temperature_present_timestep = 
    ThermodynamicTemperature::new::<degree_celsius>(143.0);

    // we won't do axial nodalisation yet...
    // but for radial nodalisation, we may perhaps use GeN-Foam code
    // now we need to relate enthalpy and temperature via some 
    // thermophysical properties
    // 
    // so all in all, we need three control volumes, 
    // one for fluid, one for shell inner node, one for shell outer 
    // node
    //
    // The shell itself will have power supplied to it and conduction 
    // heat transfer, that's all. And we want to calculate control 
    // volume enthalpy in the next timestep

    // for the outer boundary conditions, 
    // we also have an ambient_temperature
    // and an associated heat transfer coefficient


    let ambient_temperature = 
    ThermodynamicTemperature::new::<degree_celsius>(21.67);

    // let's now calculate the enthalpy at the next timestep for the 
    // outer shell layer.

    use thermal_hydraulics_rs::heat_transfer_lib::control_volume_calculations;

    // we first need timestep, and we also determine the 
    // enthalpy flows due to fluid movement to be zero
    use uom::si::time::second;
    use uom::si::power::watt;

    let timestep = Time::new::<second>(0.1);
    let solid_conductor_enthalpy_flow = 
    Power::new::<watt>(0.0);

    // the heat supplied to the system is I^2 R
    // and we know resistance is R = rho L/A
    // For the electrical heater we know potential drop across the 
    // tube is the same, therefore, V is constant
    //
    // P = V^2/R for each tube node
    // 
    // P = V^2 A_{xs} / (rho L) 
    //
    // Hence, all else equal, power scales as cross sectional area.
    // if we have heater power at 8 kW

    let total_heater_power = Power::new::<watt>(8000_f64);

    // we can take the outer node power to be the ratio of the 
    // outer node area to the whole cross sectional area
    // so circle area is pi D^2/48.0
    // 

    use uom::si::area::square_meter;

    fn circle_area(diameter: Length) -> Area {
        return diameter * diameter * PI / 4.0;
    }

    use uom::si::length::centimeter;

    // we can specify the heater inner and heater outer 
    // diameter
    let heater_od = Length::new::<centimeter>(4.0);
    let heater_id = Length::new::<centimeter>(3.81);
    let heater_midpoint = (heater_od + heater_id)/2.0;

    let heater_outer_tube_xs_area = circle_area(heater_od)
        - circle_area(heater_id);

    let heater_outer_tube_outer_node_xs_area = circle_area(heater_od)
        - circle_area(heater_midpoint);

    let heater_outer_power_fraction = heater_outer_tube_outer_node_xs_area/
        heater_outer_tube_xs_area;

    let heater_outer_power_fraction: f64 = 
    heater_outer_power_fraction.value;

    // now we can calculate heater power for outer node 

    let heater_power_outer_node: Power = 
    total_heater_power * heater_outer_power_fraction;

    // work done is zero, not considering anything

    let work_done_on_system = Power::new::<watt>(0.0);

    // actually enthalpy flow in can also be a conduction thing, 
    // but in terms of first law, it is Q to system
    // we need now to calculate heat loss to environment
    // and also heat transfer between this node and the inner node 
    //
    // Now, I'm going to assume the surface temperature is same 
    // as the finite volume temperature, though of course, one should 
    // perhaps put in a conduction resistance, so that at steady 
    // state, the solution is the same as the resistance model.
    //
    // We don't keep track of the surface temperature per se, only 
    // finite volume temperatures.
    // 
    // For heat transfer between nodes, there is also some thermal 
    // resistance between nodes 
    //

    let h_to_air: HeatTransfer = 
    HeatTransfer::new::<watt_per_square_meter_kelvin>(20.0);

    // next is the thermal conductivity
    // we'll probably need the thermal conductivity and 
    // heat capacity of steel
    //
    // From Graves,1991:
    // We write from LaTeX:

    //\begin{equation}
    //	k \left( \frac{W}{m \cdot K}  \right) = 
    //	7.9318 + 0.023051~T(K) - 6.4166*10^{-6}~T(K)^2
    //\end{equation}
    //
    //\begin{equation}
    //	c_p \left( \frac{J}{g \cdot K}  \right) = 
    //	0.4267 + 1.700* 10^{-4}~T(K) + 5.200*10^{-8}~T(K)^2
    //\end{equation}
    //
    // We may want a library or at least a function which 
    // accepts an enum saying what material it is,
    // and then based on the enum and some other things like 
    // temperature, it returns the value of the desired 
    // property in unit safe values similar to coolprop
    //
    // Though of course, we don't want to overextend ourselves
    //
    // First of course, solid properties, we'll need copper,
    // steel and fibreglass at the bare minimum

    



    let enthalpy_outer_shell_next_timestep = 
    control_volume_calculations::common_functions::
    get_control_volume_enthalpy_next_timestep(
            timestep, 
            solid_conductor_enthalpy_flow, 
            solid_conductor_enthalpy_flow, 
            heat_supplied_to_system, 
            work_done_on_system, 
            control_volume_enthalpy_current_timestep);





    


}
