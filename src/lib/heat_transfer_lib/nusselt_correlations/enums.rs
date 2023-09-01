use uom::si::{f64::*, ratio::ratio};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::input_structs::{NusseltPrandtlReynoldsData, WakaoData};

#[derive(Debug,Clone,Copy)]
pub(crate) enum NusseltCorrelation {
    PipeGnielinskiGeneric,
    PipeGnielinskiTurbulent,
    Wakao(WakaoData),
    ReynoldsPrandtl(NusseltPrandtlReynoldsData),
    PipeConstantHeatFluxFullyDeveloped,
    PipeConstantTemperatureFullyDeveloped,
}

impl NusseltCorrelation {


    /// gets the nusselt based on user choice of of correlation
    pub fn get(&self) -> Result<Ratio, ThermalHydraulicsLibError> {
        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric => todo!(),
            NusseltCorrelation::PipeGnielinskiTurbulent => todo!(),
            NusseltCorrelation::Wakao(wakao_data) => {
                return wakao_data.get();
            },
            NusseltCorrelation::ReynoldsPrandtl(reynolds_prandtl_data) => {
                return reynolds_prandtl_data.custom_reynolds_prandtl();
            },
            NusseltCorrelation::PipeConstantHeatFluxFullyDeveloped => {
                Ratio::new::<ratio>(4.354)
            },
            NusseltCorrelation::PipeConstantTemperatureFullyDeveloped => {
                Ratio::new::<ratio>(3.66)
            },
        };

        return Ok(nusselt_number);
    }

    /// Returns `true` if the nusselt correlation is [`PipeConstantHeatFlux`].
    ///
    /// [`PipeConstantHeatFlux`]: NusseltCorrelation::PipeConstantHeatFlux
    #[must_use]
    pub(crate) fn is_pipe_constant_heat_flux(&self) -> bool {
        matches!(self, Self::PipeConstantHeatFluxFullyDeveloped)
    }
}

