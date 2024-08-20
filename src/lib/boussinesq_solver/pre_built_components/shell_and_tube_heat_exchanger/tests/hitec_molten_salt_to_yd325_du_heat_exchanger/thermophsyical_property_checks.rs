
/// from 212.7 C to 279.6 C prandtl number of 
/// hitec salt should vary from 14.2 to 23.3
#[test]
pub fn prandtl_number_range_check(){

    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::{pressure::atmosphere, ratio::ratio};

    use crate::prelude::beta_testing::LiquidMaterial;

    let high_bound_temp: ThermodynamicTemperature = 
        ThermodynamicTemperature::new::<degree_celsius>(279.6);

    let low_bound_temp: ThermodynamicTemperature = 
        ThermodynamicTemperature::new::<degree_celsius>(212.7);
    let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

    let hitec = LiquidMaterial::HITEC;

    let high_bound_temp_prandtl: Ratio = 
        hitec.try_get_prandtl_liquid(high_bound_temp, atmospheric_pressure)
        .unwrap();
    let low_bound_temp_prandtl: Ratio = 
        hitec.try_get_prandtl_liquid(low_bound_temp, atmospheric_pressure)
        .unwrap();

    approx::assert_relative_eq!(
        high_bound_temp_prandtl.get::<ratio>(),
        14.2,
        max_relative=0.01
        );

    approx::assert_relative_eq!(
        low_bound_temp_prandtl.get::<ratio>(),
        23.3,
        max_relative=0.01
        );
}


/// "Plots" prandtl number from 200C to 240C in 5C intervals
#[test]
pub fn prandtl_number_range_from_200_to_240_celsius(){

    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::{pressure::atmosphere, ratio::ratio};

    use crate::prelude::beta_testing::LiquidMaterial;


    // test function
    fn check_prandtl_at_temp(temp_celsius: f64,
        prandtl_num: f64){

        let temp: ThermodynamicTemperature = 
            ThermodynamicTemperature::new::<degree_celsius>(temp_celsius);
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);
        let hitec = LiquidMaterial::HITEC;
        let calculated_prandtl: Ratio = 
            hitec.try_get_prandtl_liquid(temp, atmospheric_pressure)
            .unwrap();

            approx::assert_relative_eq!(
                calculated_prandtl.get::<ratio>(),
                prandtl_num,
                max_relative=0.01
            );

    }

    check_prandtl_at_temp(200.0, 26.24);
    check_prandtl_at_temp(205.0, 25.02);
    check_prandtl_at_temp(210.0, 23.88);
    check_prandtl_at_temp(215.0, 22.80);
    check_prandtl_at_temp(220.0, 21.78);
    check_prandtl_at_temp(225.0, 20.79);
    check_prandtl_at_temp(230.0, 19.97);
    check_prandtl_at_temp(235.0, 19.19);
    check_prandtl_at_temp(240.0, 18.46);

}


/// thermal conductivity of HITEC at 200C and 240C 
#[test]
pub fn thermal_cond_at_200_and_240_celsius(){

    use uom::si::f64::*;
    use uom::si::thermodynamic_temperature::degree_celsius;
    use uom::si::thermal_conductivity::watt_per_meter_kelvin;
    use crate::prelude::beta_testing::LiquidMaterial;


    // test function
    fn check_conductivity_at_temp(temp_celsius: f64,
        thermal_cond_watt_per_m_k: f64){

        let temp: ThermodynamicTemperature = 
            ThermodynamicTemperature::new::<degree_celsius>(temp_celsius);
        let hitec = LiquidMaterial::HITEC;
        let calculated_conductivity: ThermalConductivity = 
            hitec.try_get_thermal_conductivity(temp)
            .unwrap();

            approx::assert_relative_eq!(
                calculated_conductivity.get::<watt_per_meter_kelvin>(),
                thermal_cond_watt_per_m_k,
                max_relative=0.01
            );

    }

    check_conductivity_at_temp(200.0, 0.43602);
    check_conductivity_at_temp(240.0, 0.42806);

}
