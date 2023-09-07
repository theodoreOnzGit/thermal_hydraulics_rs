use uom::si::{f64::*, ratio::ratio};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::pipe_correlations::gnielinski_correlation_interpolated_uniform_heat_flux_liquids_developing;
/// contains information Nusselt Prandtl Reynold's
/// correlation
/// usually in the form:
///
/// Nu = a + b * Re^c * Pr^d (Pr/Pr_wall)^e
///
/// a is the constant 
/// b is the reynolds_prandtl_coefficient
/// c is the reynolds_power,
/// d is the prandtl_power,
/// e is the prandtl_correction_factor_power
#[derive(Clone,Copy,Debug, PartialEq)]
pub (in crate) struct NusseltPrandtlReynoldsData {
    pub reynolds: Ratio,
    pub prandtl_bulk: Ratio,
    pub prandtl_wall: Ratio,
    pub constant: Ratio,
    pub reynolds_prandtl_coefficient: Ratio,
    pub reynolds_power: f64,
    pub prandtl_power: f64,
    /// power for prandtl number correction factor
    pub prandtl_correction_factor_power: f64,
}

impl NusseltPrandtlReynoldsData {

    #[inline]
    pub fn custom_reynolds_prandtl(&self) 
    -> Result<Ratio,ThermalHydraulicsLibError>{
        let reynolds: Ratio =  self.reynolds;
        let prandtl_bulk: Ratio = self.prandtl_bulk;
        let prandtl_wall: Ratio = self.prandtl_wall;
        let a: Ratio = self.constant;
        let b: Ratio = self.reynolds_prandtl_coefficient;
        let c: f64 = self.reynolds_power;
        let d: f64 = self.prandtl_power;
        let e: f64 = self.prandtl_correction_factor_power;

        let nusselt: Ratio = a + 
        b * reynolds.get::<ratio>().powf(c) 
        * prandtl_bulk.get::<ratio>().powf(d) 
        * (prandtl_bulk/prandtl_wall).get::<ratio>().powf(e);

        return Ok(nusselt);
    }

    #[inline]
    pub fn ciet_version_2_heater_uncorrected(&self) -> 
    Result<Ratio, ThermalHydraulicsLibError>{

        let reynolds = self.reynolds;
        let prandtl = self.prandtl_bulk;
        let reynolds_power_0_836 = reynolds.value.powf(0.836);
        let prandtl_power_0_333 = prandtl.value.powf(0.333333333333333);

        let nusselt = Ratio::new::<ratio>(
            0.04179 * reynolds_power_0_836 * prandtl_power_0_333);

        Ok(nusselt)
    }

    #[inline]
    pub fn ciet_version_2_heater_prandtl_corrected(&self) -> 
    Result<Ratio, ThermalHydraulicsLibError>{
        let nusself_uncorrected 
        =  {
            let ref this = self;

            let reynolds = this.reynolds;
            let prandtl = this.prandtl_bulk;
            let reynolds_power_0_836 = reynolds.value.powf(0.836);
            let prandtl_power_0_333 = prandtl.value.powf(0.333333333333333);

            let nusselt_uncorrected = Ratio::new::<ratio>(
                0.04179 * reynolds_power_0_836 * prandtl_power_0_333);

            nusselt_uncorrected
        };

        let prandtl_wall = self.prandtl_wall;
        let prandtl_bulk = self.prandtl_bulk;

        let prandtl_bulk_to_wall_ratio = prandtl_wall/prandtl_bulk;

        let correction_factor: f64 
        = prandtl_bulk_to_wall_ratio.get::<ratio>().powf(0.11);

        return Ok(nusself_uncorrected*correction_factor);

    }
}


/// Wakao, N., & Funazkri, T. (1978). Effect 
/// of fluid dispersion coefficients on particle-to-fluid mass 
/// transfer coefficients in packed beds: correlation of 
/// Sherwood numbers. Chemical Engineering Science, 33(10), 1375-1384.
#[derive(Clone,Copy,Debug, PartialEq)]
pub (in crate) struct WakaoData {
    pub reynolds: Ratio,
    pub prandtl_bulk: Ratio,
}

impl WakaoData {

    #[inline]
    pub fn get(&self) 
    -> Result<Ratio,ThermalHydraulicsLibError>{
        let reynolds: Ratio =  self.reynolds;
        let prandtl_bulk: Ratio = self.prandtl_bulk;
        let a: Ratio = Ratio::new::<ratio>(2.0);
        let b: Ratio = Ratio::new::<ratio>(1.1);
        let c: f64 = 0.3333333333;
        let d: f64 = 0.6;

        let nusselt: Ratio = a + 
        b * reynolds.get::<ratio>().powf(c) 
        * prandtl_bulk.get::<ratio>().powf(d);

        return Ok(nusselt);
    }
}


/// contains data for gnielinski 
/// correlation of various
#[derive(Clone,Copy,Debug, PartialEq)]
pub (in crate) struct GnielinskiData {
    pub reynolds: Ratio,
    pub prandtl_bulk: Ratio,
    pub prandtl_wall: Ratio,
    pub darcy_friction_factor: Ratio,
    pub length_to_diameter: Ratio
}

impl GnielinskiData {

    #[inline]
    pub fn get_nusselt_for_developing_flow(&self) 
    -> Result<Ratio,ThermalHydraulicsLibError>{
        let reynolds: Ratio =  self.reynolds;
        let prandtl_bulk: Ratio = self.prandtl_bulk;
        let prandtl_wall: Ratio = self.prandtl_wall;
        let darcy_friction_factor = self.darcy_friction_factor;
        let length_to_diameter = self.length_to_diameter;

        let nusselt_value = 
        gnielinski_correlation_interpolated_uniform_heat_flux_liquids_developing(
            reynolds.get::<ratio>(),
            prandtl_bulk.get::<ratio>(),
            prandtl_wall.get::<ratio>(),
            darcy_friction_factor.get::<ratio>(),
            length_to_diameter.get::<ratio>(),
        );

        return Ok(
            Ratio::new::<ratio>(nusselt_value)
        );

    }

}

