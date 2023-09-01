use uom::si::{f64::*, ratio::ratio};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
/// contains information Nusselt Prandtl Reynold's
/// correlation
/// usually in the form:
///
/// Nu = a + b * Re^c * Pr^d (Pr/Pr_wall)^e
#[derive(Clone,Copy,Debug, PartialEq)]
pub (in crate) struct NusseltPrandtlReynoldsData {
    pub reynolds: Ratio,
    pub prandtl_bulk: Ratio,
    pub prandtl_wall: Ratio,
    pub a: Ratio,
    pub b: Ratio,
    pub c: f64,
    pub d: f64,
    pub e: f64,
}

impl NusseltPrandtlReynoldsData {

    #[inline]
    pub fn custom_reynolds_prandtl(&self) 
    -> Result<Ratio,ThermalHydraulicsLibError>{
        let reynolds: Ratio =  self.reynolds;
        let prandtl_bulk: Ratio = self.prandtl_bulk;
        let prandtl_wall: Ratio = self.prandtl_wall;
        let a: Ratio = self.a;
        let b: Ratio = self.b;
        let c: f64 = self.c;
        let d: f64 = self.d;
        let e: f64 = self.e;

        let nusselt: Ratio = a + 
        b * reynolds.get::<ratio>().powf(c) 
        * prandtl_bulk.get::<ratio>().powf(d) 
        * (prandtl_bulk/prandtl_wall).get::<ratio>().powf(e);

        return Ok(nusselt);
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


#[derive(Clone,Copy,Debug, PartialEq)]
pub (in crate) struct Gnielinski {
    pub reynolds: Ratio,
    pub prandtl_bulk: Ratio,
    pub prandtl_wall: Ratio,
    pub a: Ratio,
    pub b: Ratio,
    pub c: f64,
    pub d: f64,
    pub e: f64,
    pub friction_factor: Ratio,
    pub length_to_diameter: Ratio
}

impl Gnielinski {

    #[inline]
    pub fn custom_reynolds_prandtl(&self) 
    -> Result<Ratio,ThermalHydraulicsLibError>{
        let reynolds: Ratio =  self.reynolds;
        let prandtl_bulk: Ratio = self.prandtl_bulk;
        let prandtl_wall: Ratio = self.prandtl_wall;
        let a: Ratio = self.a;
        let b: Ratio = self.b;
        let c: f64 = self.c;
        let d: f64 = self.d;
        let e: f64 = self.e;

        let nusselt: Ratio = a + 
        b * reynolds.get::<ratio>().powf(c) 
        * prandtl_bulk.get::<ratio>().powf(d) 
        * (prandtl_bulk/prandtl_wall).get::<ratio>().powf(e);

        return Ok(nusselt);
    }
}
