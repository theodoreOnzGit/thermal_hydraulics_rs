use uom::si::f64::*;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;

// contains specific enthalpy data for all materials
mod enthalpy_data;
use enthalpy_data::*;

// contains functions to map enthalpy to a temperature 
mod temperature_from_specific_enthalpy;
use temperature_from_specific_enthalpy::*;

/// returns specific enthaply for a given material 
/// specific_enthalpy is defined as 0 for 0 degree_celsius
/// for any material, that is 273.15 K
///
/// ```rust 
/// use uom::si::f64::*;
/// use uom::si::specific_heat_capacity::{joule_per_kilogram_kelvin,
/// joule_per_gram_degree_celsius};
/// use uom::si::thermodynamic_temperature::kelvin;
/// use uom::si::temperature_interval::degree_celsius;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// SolidMaterial::{SteelSS304L,Copper};
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// Material;
/// use thermal_hydraulics_rs::heat_transfer_lib::thermophysical_properties::
/// specific_enthalpy::specific_enthalpy;
///
/// use uom::si::pressure::atmosphere;
///
/// let steel = Material::Solid(SteelSS304L);
/// let steel_temp = ThermodynamicTemperature::new::<kelvin>(273.15);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// // enthalpy should be zero at 273.15 K
///
/// let steel_enthalpy_273_15_kelvin = 
/// specific_enthalpy(steel, steel_temp, pressure);
///
/// approx::assert_relative_eq!(
///     0.0,
///     steel_enthalpy_273_15_kelvin.unwrap().value,
///     max_relative=0.045);
/// 
/// // we can also calculate enthalpy change of copper 
/// // from 375K to 425K
/// let test_temperature_1 = ThermodynamicTemperature::new::
/// <kelvin>(375.0);
/// let test_temperature_2 = ThermodynamicTemperature::new::
/// <kelvin>(425.0);
///
/// let copper = Material::Solid(Copper);
///
/// let copper_enthalpy_change = 
/// specific_enthalpy(copper, test_temperature_2, pressure).unwrap()
/// - specific_enthalpy(copper, test_temperature_1, pressure).unwrap();
///
/// // http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/sphtt.html
/// // https://www.engineeringtoolbox.com/specific-heat-metals-d_152.html
/// // copper at 20C has heat capacity of 
/// // 0.386 J/(g K)
/// // going to use this to estimate a ballpark figure to find enthalpy 
/// // h = cp(T2 - T1)
/// 
/// // we can't usually subtract thermodynamic temperatures from each 
/// // other, we need a termpature interval
/// // 
///
/// let cp_copper_20_c = 
/// SpecificHeatCapacity::new::<joule_per_gram_degree_celsius>(0.386);
/// 
/// let temperature_difference = 
/// TemperatureInterval::new::<degree_celsius>(
/// test_temperature_2.value - test_temperature_1.value);
///
/// let specific_enthalpy_ballpark = 
/// cp_copper_20_c * temperature_difference;
/// 
/// // the ballpark value is 19300 J/kg
/// approx::assert_relative_eq!(
///     specific_enthalpy_ballpark.value,
///     19300.0,
///     max_relative=0.0001);
///
/// // it's less than 4% different from the ballpark value
/// // This means the copper enthalpy change should be quite reasonable
///
/// approx::assert_relative_eq!(
///     specific_enthalpy_ballpark.value,
///     copper_enthalpy_change.value,
///     max_relative=0.04);
///
/// ``` 
pub fn specific_enthalpy(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<AvailableEnergy, String> {

    let specific_enthalpy: AvailableEnergy = match material {
        Material::Solid(_) => solid_specific_enthalpy(material, temperature),
        Material::Liquid(_) => liquid_specific_enthalpy(material, temperature)
    };

    return Ok(specific_enthalpy);
}




