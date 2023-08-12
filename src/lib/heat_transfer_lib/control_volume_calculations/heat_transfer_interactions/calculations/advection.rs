use uom::si::f64::*;
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

