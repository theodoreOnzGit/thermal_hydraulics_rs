use uom::si::{f64::*, ratio::ratio};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::input_structs::{NusseltPrandtlReynoldsData, WakaoData, GnielinskiData};

#[derive(Debug,Clone,Copy,Default, PartialEq)]
pub(crate) enum NusseltCorrelation {
    PipeGnielinskiGeneric(GnielinskiData),
    PipeGnielinskiTurbulent(GnielinskiData),
    Wakao(WakaoData),
    ReynoldsPrandtl(NusseltPrandtlReynoldsData),
    #[default]
    PipeConstantHeatFluxFullyDeveloped,
    PipeConstantTemperatureFullyDeveloped,
    CIETHeaterVersion2(NusseltPrandtlReynoldsData),
}

impl NusseltCorrelation {


    /// gets the nusselt based on user choice of of correlation
    pub fn get(&self) -> Result<Ratio, ThermalHydraulicsLibError> {
        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric(data) => {
                return data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::PipeGnielinskiTurbulent(data) => {
                return data.get_nusselt_for_developing_flow();
            },
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
            NusseltCorrelation::CIETHeaterVersion2(data) => {
                return data.ciet_version_2_heater_prandtl_corrected();
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

