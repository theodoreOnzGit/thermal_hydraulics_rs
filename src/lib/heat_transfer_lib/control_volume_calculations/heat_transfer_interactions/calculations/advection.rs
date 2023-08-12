use uom::si::{f64::*, mass_rate::kilogram_per_second, available_energy::joule_per_kilogram};
/// now, advection is quite tricky because for conduction, the 
/// heat transfer formula is for two control volumes cv_a and cv_b
/// can be as follows 
///
/// (cv_a) --------------- (cv_b)
///  
///  T_a                    T_b 
///
///  Q_(ab) = -H(T_b - T_a)
/// 
/// Q_(ab) is heat transfer rate (watts) from a to b
/// H is conductance, not heat transfer coefficient
/// it has units of watts kelvin
///
/// For advection in contrast, it depends on flow
/// 
/// (cv_a) --------------- (cv_b)
///  
///  T_a                    T_b 
///
/// For flow from a to b:
/// Q_(ab) = m h(T_a)
///
/// For flow from b to a 
/// Q_(ab) = -m h(T_b)
///
/// Here, the enthalpy transfer only depends on one of the body's 
/// temperature, which is directly dependent on mass flow 
#[inline]
pub fn advection_heat_rate(mass_flow_from_a_to_b: MassRate,
    specific_enthalpy_of_a: AvailableEnergy,
    specific_enthalpy_of_b: AvailableEnergy,) -> Result<Power, String> {

    // if mass flow from a to b is less than 0 
    // then mass flows from b to a
    //
    // power is dependent on control volume b
    let heat_rate: Power;
    if mass_flow_from_a_to_b.value < 0.0_f64 {
        heat_rate = mass_flow_from_a_to_b * specific_enthalpy_of_b;
    } else {
        heat_rate = mass_flow_from_a_to_b * specific_enthalpy_of_a;
    }

    Ok(heat_rate)
}

/// this is just a test for me to check if the advection signs work 
/// correctly
///
/// Suppose there is this following scenario:
///
/// (cv_a) --------------- (cv_b)
///  
///  h_a                    h_b 
///
///  mass flow from a to b is 5 kg/s 
///
///  h_a = 350 J/kg
///  h_b = 550 J/kg
///
///
///
#[test]
pub fn test_advection_signs(){

    let forward_mass_flow: MassRate = 
    MassRate::new::<kilogram_per_second>(5.0);
    let backwards_mass_flow: MassRate = 
    MassRate::new::<kilogram_per_second>(-5.0);

    let enthalpy_a: AvailableEnergy = 
    AvailableEnergy::new::<joule_per_kilogram>(350.0);
    let enthalpy_b: AvailableEnergy = 
    AvailableEnergy::new::<joule_per_kilogram>(550.0);

    // if forward flow from a to b, then expected power in forward 
    // flow is
    let expected_forward_flow_enthalpy_a_to_b: Power 
    = forward_mass_flow * enthalpy_a;

    let test_forward_flow_enthalpy_rate: Power 
    = advection_heat_rate(forward_mass_flow,
        enthalpy_a,
        enthalpy_b).unwrap();

    assert_eq!(test_forward_flow_enthalpy_rate, 
        expected_forward_flow_enthalpy_a_to_b);

    // if backward flow from b to a, expected power in backward flow is:
    let expected_forward_flow_enthalpy_a_to_b: Power 
    = backwards_mass_flow * enthalpy_b;

    let test_backward_flow_enthalpy_rate: Power 
    = advection_heat_rate(backwards_mass_flow,
        enthalpy_a,
        enthalpy_b).unwrap();

    // check if flows are different, for forward and backward 
    // and check if flows are correct
    // 
    assert_ne!(test_forward_flow_enthalpy_rate, 
        expected_forward_flow_enthalpy_a_to_b);
    assert_eq!(test_backward_flow_enthalpy_rate, 
        expected_forward_flow_enthalpy_a_to_b);

}
