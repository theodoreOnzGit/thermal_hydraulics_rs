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
