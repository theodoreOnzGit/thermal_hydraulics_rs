use uom::si::f64::MassDensity;
use uom::si::f64::Pressure;
use uom::si::f64::ThermodynamicTemperature;
use uom::si::mass_density::kilogram_per_cubic_meter;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::liquid_database;
use super::liquid_database::flibe::get_flibe_density;
use super::liquid_database::hitec_nitrate_salt::get_hitec_density;
use super::liquid_database::yd_325_heat_transfer_oil::get_yd325_density;
use super::solid_database::custom_solid_material;
use super::LiquidMaterial;
use super::Material;
use super::SolidMaterial;
use super::SolidMaterial::*;
use super::LiquidMaterial::*;
use super::liquid_database::dowtherm_a::get_dowtherm_a_density;

/// returns a density given a material, temperature and pressure
///
/// example:
///
/// ```rust
/// use uom::si::f64::*;
/// use uom::si::pressure::atmosphere;
/// use uom::si::thermodynamic_temperature::kelvin;
/// use thermal_hydraulics_rs::boussinesq_solver::
/// boussinesq_thermophysical_properties::density::try_get_rho;
///
/// use thermal_hydraulics_rs::boussinesq_solver::
/// boussinesq_thermophysical_properties::SolidMaterial::SteelSS304L;
///
/// use thermal_hydraulics_rs::boussinesq_solver::
/// boussinesq_thermophysical_properties::Material;
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
    _pressure: Pressure) -> Result<MassDensity, ThermalHydraulicsLibError> {

    let density: MassDensity = match material {
        Material::Solid(_) => solid_density(material, temperature)?,
        Material::Liquid(_) => liquid_density(material, temperature)?
    };

    return Ok(density);
}

impl Material {
    /// returns density of the material
    pub fn density(&self,
        temperature: ThermodynamicTemperature,
        _pressure: Pressure) -> Result<MassDensity, ThermalHydraulicsLibError>{

    let density: MassDensity = match self {
        Material::Solid(_) => solid_density(self.clone(), temperature)?,
        Material::Liquid(_) => liquid_density(self.clone(), temperature)?
    };

    return Ok(density);
        

    }

}

// should the material happen to be a solid, use this function
fn solid_density(material: Material,
    solid_temp: ThermodynamicTemperature) -> Result<MassDensity,ThermalHydraulicsLibError>{

    // first match the enum

    let solid_material: SolidMaterial = match material {
        Material::Solid(SteelSS304L) => SteelSS304L,
        Material::Solid(Fiberglass) => Fiberglass,
        Material::Solid(Copper) => Copper,
        Material::Solid(CustomSolid((low_bound_temp,high_bound_temp),cp,k,rho,roughness)) => {
            CustomSolid((low_bound_temp,high_bound_temp), cp, k, rho,roughness)
        },
        Material::Liquid(_) => {
            println!("solid_density, use SolidMaterial enums only");
            return Err(ThermalHydraulicsLibError::TypeConversionErrorMaterial);
        }
    };

    let density: MassDensity = match solid_material {
        Fiberglass => fiberglass_density()?,
        SteelSS304L => steel_ss_304_l_density()?,
        Copper => copper_density()?,
        CustomSolid((low_bound_temp,high_bound_temp),_cp,_k,rho_fn,_roughness) => {
            custom_solid_material::get_custom_solid_density(
                solid_temp, 
                rho_fn, 
                high_bound_temp, 
                low_bound_temp)?
        },
    };

    return Ok(density);


}




// should the material happen to be a liquid, use this function
fn liquid_density(material: Material, 
    fluid_temp: ThermodynamicTemperature) -> Result<MassDensity,ThermalHydraulicsLibError> {

    let liquid_material: LiquidMaterial = match material {
        Material::Liquid(DowthermA) => DowthermA,
        Material::Liquid(TherminolVP1) => TherminolVP1,
        Material::Liquid(HITEC) => HITEC,
        Material::Liquid(YD325) => YD325,
        Material::Liquid(FLiBe) => FLiBe,
        Material::Liquid(CustomLiquid((low_bound_temp,high_bound_temp),cp,k,mu,rho)) => {
            CustomLiquid((low_bound_temp,high_bound_temp), cp, k, mu, rho)
        },

        Material::Solid(_) => panic!("liquid_density, use LiquidMaterial enums only")
    };

    let density: MassDensity = match liquid_material {
        DowthermA => dowtherm_a_density(fluid_temp)?,
        TherminolVP1 => dowtherm_a_density(fluid_temp)?,
        HITEC => get_hitec_density(fluid_temp)?,
        YD325 => get_yd325_density(fluid_temp)?,
        FLiBe => get_flibe_density(fluid_temp)?,
        CustomLiquid((low_bound_temp,high_bound_temp), _cp, _k, _mu, rho_fn) => {
            liquid_database::custom_liquid_material
                ::get_custom_fluid_density(fluid_temp, 
                    rho_fn, 
                    high_bound_temp, 
                    low_bound_temp)?
        },
    };

    return Ok(density);
}

impl LiquidMaterial {

    /// returns density of liquid material
    pub fn density(&self,
        fluid_temp: ThermodynamicTemperature,) -> 
    Result<MassDensity,ThermalHydraulicsLibError> {

        let density: MassDensity = match &self.clone() {
            DowthermA => dowtherm_a_density(fluid_temp)?,
            TherminolVP1 => dowtherm_a_density(fluid_temp)?,
            HITEC => get_hitec_density(fluid_temp)?,
            YD325 => get_yd325_density(fluid_temp)?,
            FLiBe => get_flibe_density(fluid_temp)?,
            CustomLiquid((low_bound_temp,high_bound_temp), _cp, _k, _mu, rho_fn) => {
                liquid_database::custom_liquid_material
                    ::get_custom_fluid_density(fluid_temp, 
                        *rho_fn, 
                        *high_bound_temp, 
                        *low_bound_temp)?
            },
        };

        Ok(density)

    }
}

#[inline]
fn fiberglass_density() -> Result<MassDensity,ThermalHydraulicsLibError> {
    // density ranges not quite given in original text 
    // Zou, Ling, Rui Hu, and Anne Charpentier. SAM code validation 
    // using the compact integral effects test (CIET) experimental data. 
    // No. ANL/NSE-19/11. 
    // Argonne National Lab.(ANL), Argonne, IL (United States), 2019.
    return Ok(MassDensity::new::<kilogram_per_cubic_meter>(20.0));
}

#[inline]
fn steel_ss_304_l_density() -> Result<MassDensity,ThermalHydraulicsLibError> {
    // density ranges not quite given in original text 
    // Zou, Ling, Rui Hu, and Anne Charpentier. SAM code validation 
    // using the compact integral effects test (CIET) experimental data. 
    // No. ANL/NSE-19/11. 
    // Argonne National Lab.(ANL), Argonne, IL (United States), 2019.
    return Ok(MassDensity::new::<kilogram_per_cubic_meter>(8030.0));
}

#[inline]
fn copper_density() -> Result<MassDensity,ThermalHydraulicsLibError> {
    // density ranges not quite given in original text 
    // Zou, Ling, Rui Hu, and Anne Charpentier. SAM code validation 
    // using the compact integral effects test (CIET) experimental data. 
    // No. ANL/NSE-19/11. 
    // Argonne National Lab.(ANL), Argonne, IL (United States), 2019.
    return Ok(MassDensity::new::<kilogram_per_cubic_meter>(8940.0));
}

#[inline]
fn dowtherm_a_density(fluid_temp: ThermodynamicTemperature) -> 
Result<MassDensity,ThermalHydraulicsLibError>{
    return get_dowtherm_a_density(fluid_temp);
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
