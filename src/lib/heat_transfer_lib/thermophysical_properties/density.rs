use uom::si::f64::MassDensity;
use uom::si::f64::Pressure;
use uom::si::f64::ThermodynamicTemperature;
use uom::si::mass_density::kilogram_per_cubic_meter;
use crate::fluid_mechanics_lib::therminol_component::
dowtherm_a_properties::getDowthermADensity;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;

/// returns a density given a material, temperature and pressure
///
/// example:
///
/// ```rust
/// use uom::si::f64::*;
/// use uom::si::pressure::atmosphere;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::density::try_get_rho;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::SolidMaterial::SteelSS304L;
///
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// thermophysical_properties::Material;
///
/// let steel = Material::Solid(SteelSS304L);
/// let temperature = ThermodynamicTemperature::new::<kelvin>(396.0);
/// let pressure = Pressure::new::<atmosphere>(1.0);
///
/// let density_result = try_get_rho(steel, temperature, pressure);
///
/// 
/// ```
#[inline]
pub fn try_get_rho(material: Material, 
    temperature: ThermodynamicTemperature,
    _pressure: Pressure) -> Result<MassDensity, String> {

    let density: MassDensity = match material {
        Material::Solid(_) => solid_density(material, temperature),
        Material::Liquid(_) => liquid_density(material, temperature)
    };

    return Ok(density);
}

impl Material {
    /// returns density of the material
    pub fn density(&self,
        temperature: ThermodynamicTemperature,
        _pressure: Pressure) -> Result<MassDensity, ThermalHydraulicsLibError>{

    let density: MassDensity = match self {
        Material::Solid(_) => solid_density(self.clone(), temperature),
        Material::Liquid(_) => liquid_density(self.clone(), temperature)
    };

    return Ok(density);
        

    }

}

// should the material happen to be a solid, use this function
fn solid_density(material: Material,
    _temperature: ThermodynamicTemperature) -> MassDensity{

    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Liquid(_) => panic!("solid_density, use SolidMaterial enums only")
    };

    let density: MassDensity = match solid_material {
        Fiberglass => fiberglass_density() ,
        SteelSS304L => steel_ss_304_l_density(),
        Copper => copper_density(),
    };

    return density;


}




// should the material happen to be a liquid, use this function
fn liquid_density(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> MassDensity {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Solid(_) => panic!("liquid_density, use LiquidMaterial enums only")
    };

    let density: MassDensity = match liquid_material {
        DowthermA => dowtherm_a_density(fluid_temp),
        TherminolVP1 => dowtherm_a_density(fluid_temp)
    };

    return density;
}

impl LiquidMaterial {

    /// returns density of liquid material
    pub fn density(&self,
        fluid_temp: ThermodynamicTemperature,) -> 
    Result<MassDensity,ThermalHydraulicsLibError> {

        let density: MassDensity = match &self.clone() {
            DowthermA => dowtherm_a_density(fluid_temp),
            TherminolVP1 => dowtherm_a_density(fluid_temp)
        };

        Ok(density)

    }
}

#[inline]
fn fiberglass_density() -> MassDensity {
    return MassDensity::new::<kilogram_per_cubic_meter>(20.0);
}

#[inline]
fn steel_ss_304_l_density() -> MassDensity {
    return MassDensity::new::<kilogram_per_cubic_meter>(8030.0);
}

#[inline]
fn copper_density() -> MassDensity {
    return MassDensity::new::<kilogram_per_cubic_meter>(8940.0);
}

#[inline]
fn dowtherm_a_density(fluid_temp: ThermodynamicTemperature) -> MassDensity{
    return getDowthermADensity(fluid_temp);
}

#[test]
pub fn density_test_steel(){

    use uom::si::thermodynamic_temperature::kelvin;
    use uom::si::pressure::atmosphere;
    let steel = Material::Solid(SteelSS304L);
    let temperature = ThermodynamicTemperature::new::<kelvin>(396.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let density = try_get_rho(steel, temperature, pressure);

    approx::assert_relative_eq!(
        8030_f64,
        density.unwrap().value,
        max_relative=0.01);
}
