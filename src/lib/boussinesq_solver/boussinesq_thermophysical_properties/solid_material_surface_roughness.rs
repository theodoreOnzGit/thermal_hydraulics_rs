use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use uom::si::{f64::*, length::{micrometer, millimeter}};

use super::{SolidMaterial, Material};

impl SolidMaterial {

    /// returns surface roughness for various materials
    pub fn surface_roughness(&self) -> Result<Length,ThermalHydraulicsLibError> {
        let roughness: Length = match self {
            // Value from: Perry's chemical Engineering handbook 
            // 8th edition Table 6-1 
            // commercial steel or wrought iron 
            // Perry, R. H., & DW, G. (2007). 
            // Perry’s chemical engineers’ handbook, 
            // 8th illustrated ed. New York: McGraw-Hill.
            SolidMaterial::SteelSS304L => {
                Length::new::<millimeter>(0.0457)
            },
            // Arenales, M. R. M., Kumar, S., 
            // Kuo, L. S., & Chen, P. H. (2020). 
            // Surface roughness variation effects on copper tubes in 
            // pool boiling of water. International Journal of 
            // Heat and Mass Transfer, 151, 119399.
            SolidMaterial::Copper => {
                Length::new::<micrometer>(0.544)
            },
            // Value from: Perry's chemical Engineering handbook 
            // 8th edition Table 6-1 
            // generic value for drawn tubing
            // Perry, R. H., & DW, G. (2007). 
            // Perry’s chemical engineers’ handbook, 
            // 8th illustrated ed. New York: McGraw-Hill.
            SolidMaterial::Fiberglass => {
                Length::new::<millimeter>(0.00152)
            },
        };

        Ok(roughness)
    }
}

impl Material {
    /// wrapper to help return surface roughness
    pub fn surface_roughness(&self) -> Result<Length,ThermalHydraulicsLibError>{
        match self {
            Material::Solid(solid_material) => {
                return solid_material.surface_roughness();
            },
            Material::Liquid(_) => {
                Err(ThermalHydraulicsLibError::TypeConversionErrorMaterial)
            },
        }
    }
}
