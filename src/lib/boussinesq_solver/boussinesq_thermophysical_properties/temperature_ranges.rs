use peroxide::fuga::P;
use uom::si::f64::*;

use super::{LiquidMaterial, Material, SolidMaterial};

impl Material {

    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn max_temperature(&self) -> ThermodynamicTemperature {
        match self {
            Material::Solid(solid) => {
                solid.max_temperature()
            },
            Material::Liquid(liquid) => {
                liquid.max_temperature()
            },
        }
    }
    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn min_temperature(&self) -> ThermodynamicTemperature {

        match self {
            Material::Solid(solid) => {
                solid.min_temperature()
            },
            Material::Liquid(liquid) => {
                liquid.min_temperature()
            },
        }
    }
}

impl LiquidMaterial {

    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn max_temperature(&self) -> ThermodynamicTemperature {
        match self {
            LiquidMaterial::TherminolVP1 => todo!(),
            LiquidMaterial::DowthermA => todo!(),
            LiquidMaterial::HITEC => todo!(),
            LiquidMaterial::YD325 => todo!(),
            LiquidMaterial::FLiBe => todo!(),
            LiquidMaterial::FLiNaK => todo!(),
            LiquidMaterial::CustomLiquid((_lower_bound, upper_bound)
                , _, _, _, _) => {
                *upper_bound
            },
        }
    }
    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn min_temperature(&self) -> ThermodynamicTemperature {
        match self {
            LiquidMaterial::TherminolVP1 => todo!(),
            LiquidMaterial::DowthermA => todo!(),
            LiquidMaterial::HITEC => todo!(),
            LiquidMaterial::YD325 => todo!(),
            LiquidMaterial::FLiBe => todo!(),
            LiquidMaterial::FLiNaK => todo!(),
            LiquidMaterial::CustomLiquid((lower_bound, _upper_bound)
                , _, _, _, _) => {
                *lower_bound
            },
        }

    }
}
impl SolidMaterial {

    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn max_temperature(&self) -> ThermodynamicTemperature {
        match self {
            SolidMaterial::SteelSS304L => todo!(),
            SolidMaterial::Copper => todo!(),
            SolidMaterial::Fiberglass => todo!(),
            SolidMaterial::CustomSolid((_lower_bound,upper_bound), 
                _, _, _, _) => {
                *upper_bound
            },
        }
    }
    /// gives the maximum temperature for the correlations in the
    /// given material
    pub fn min_temperature(&self) -> ThermodynamicTemperature {
        match self {
            SolidMaterial::SteelSS304L => todo!(),
            SolidMaterial::Copper => todo!(),
            SolidMaterial::Fiberglass => todo!(),
            SolidMaterial::CustomSolid((lower_bound, _upper_bound), 
                _, _, _, _) => {
                *lower_bound
            },
        }

    }
}
