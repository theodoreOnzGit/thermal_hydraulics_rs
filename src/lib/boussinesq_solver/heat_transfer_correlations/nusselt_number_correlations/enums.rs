use uom::si::{f64::*, ratio::ratio};

use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use super::input_structs::{NusseltPrandtlReynoldsData, WakaoData, GnielinskiData};

/// Contains a collection of nusselt number correlations for use 
///
/// still under experimentation, so quite unstable
///
/// for now, you are supposed to construct a struct containing Pr and Re 
/// and etc and fit them into the enum,
///
/// use the get method to obtain the Nusselt Number
#[derive(Debug,Clone,Copy,Default, PartialEq)]
pub enum NusseltCorrelation {

    /// pipe nusselt number using Gnielinski Correlation 
    /// for laminar, turbulent and transition region
    ///
    /// laminar flow assumes constant heat flux
    PipeGnielinskiGeneric(GnielinskiData),

    /// nusselt number only for turbulent
    /// flow in pipes
    PipeGnielinskiTurbulent(GnielinskiData),

    /// nusselt number for porous media 
    /// especially packed beds
    /// based on Wakao Correlation
    ///
    /// note: reynolds number based on pebble diameter
    /// Wakao, N., & Funazkri, T. (1978). Effect 
    /// of fluid dispersion coefficients on particle-to-fluid mass 
    /// transfer coefficients in packed beds: correlation of 
    /// Sherwood numbers. Chemical Engineering Science, 33(10), 1375-1384.
    Wakao(WakaoData),

    /// generic reynolds prandtl power correlation
    /// usually in the form:
    ///
    /// Nu = a + b * Re^c * Pr^d (Pr/Pr_wall)^e
    ///
    /// a is the constant 
    /// b is the reynolds_prandtl_coefficient
    /// c is the reynolds_power,
    /// d is the prandtl_power,
    /// e is the prandtl_correction_factor_power
    ReynoldsPrandtl(NusseltPrandtlReynoldsData),

    /// returns a nusselt number of 4.36 for fully developed 
    /// constant heat flux flow
    #[default]
    PipeConstantHeatFluxFullyDeveloped,
    /// returns a nusselt number of 3.66 for fully developed 
    /// constant temperature flow
    PipeConstantTemperatureFullyDeveloped,

    /// ciet heater correlation for version 2, 
    ///
    /// Nu = 0.04179 * reynolds^0.836 * prandtl^0.333
    /// * (Pr_wall/Pr_bulk)^0.11
    CIETHeaterVersion2(NusseltPrandtlReynoldsData),
}

impl NusseltCorrelation {


    /// gets the nusselt based on user choice of of correlation
    pub fn try_get(&self) -> Result<Ratio, ThermalHydraulicsLibError> {
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

    /// gets an estimate for the nusselt number based on user choice 
    /// of correlation, ignores wall temperature 
    pub fn estimate_based_on_prandtl_and_reynolds(&self,
    prandtl_number_input: Ratio,
    reynolds_number_input: Ratio,) -> Result<Ratio, ThermalHydraulicsLibError>{

        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric(data) => {

                let mut modified_data = data.clone();
                modified_data.prandtl_wall = prandtl_number_input;
                modified_data.prandtl_bulk = prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::PipeGnielinskiTurbulent(data) => {
                let mut modified_data = data.clone();
                modified_data.prandtl_wall = prandtl_number_input;
                modified_data.prandtl_bulk = prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::Wakao(wakao_data) => {
                let mut modified_data = wakao_data.clone();
                modified_data.prandtl_bulk = prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get();
            },
            NusseltCorrelation::ReynoldsPrandtl(reynolds_prandtl_data) => {
                let mut modified_data = reynolds_prandtl_data.clone();
                modified_data.prandtl_wall = prandtl_number_input;
                modified_data.prandtl_bulk = prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;

                return modified_data.custom_reynolds_prandtl();
            },
            NusseltCorrelation::PipeConstantHeatFluxFullyDeveloped => {
                Ratio::new::<ratio>(4.354)
            },
            NusseltCorrelation::PipeConstantTemperatureFullyDeveloped => {
                Ratio::new::<ratio>(3.66)
            },
            NusseltCorrelation::CIETHeaterVersion2(data) => {
                let mut modified_data = data.clone();
                modified_data.prandtl_bulk = prandtl_number_input;
                modified_data.prandtl_wall = prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;

                return modified_data.ciet_version_2_heater_prandtl_corrected();
            },
        };

        return Ok(nusselt_number);
    }

    /// gets an estimate for the nusselt number based on user choice 
    /// of correlation, ignores wall temperature 
    pub fn estimate_based_on_prandtl_reynolds_and_wall_correction(&self,
    bulk_prandtl_number_input: Ratio,
    wall_prandtl_number_input: Ratio,
    reynolds_number_input: Ratio,) -> Result<Ratio, ThermalHydraulicsLibError>{

        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric(data) => {

                let mut modified_data = data.clone();
                modified_data.prandtl_wall = wall_prandtl_number_input;
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::PipeGnielinskiTurbulent(data) => {
                let mut modified_data = data.clone();
                modified_data.prandtl_wall = wall_prandtl_number_input;
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::Wakao(wakao_data) => {
                let mut modified_data = wakao_data.clone();
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get();
            },
            NusseltCorrelation::ReynoldsPrandtl(reynolds_prandtl_data) => {
                let mut modified_data = reynolds_prandtl_data.clone();
                modified_data.prandtl_wall = wall_prandtl_number_input;
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;

                return modified_data.custom_reynolds_prandtl();
            },
            NusseltCorrelation::PipeConstantHeatFluxFullyDeveloped => {
                Ratio::new::<ratio>(4.354)
            },
            NusseltCorrelation::PipeConstantTemperatureFullyDeveloped => {
                Ratio::new::<ratio>(3.66)
            },
            NusseltCorrelation::CIETHeaterVersion2(data) => {
                let mut modified_data = data.clone();
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.prandtl_wall = wall_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;

                return modified_data.ciet_version_2_heater_prandtl_corrected();
            },
        };

        return Ok(nusselt_number);
    }

    /// Returns `true` if the nusselt correlation is [`PipeConstantHeatFlux`].
    ///
    /// [`PipeConstantHeatFlux`]: NusseltCorrelation::PipeConstantHeatFlux
    #[must_use]
    pub(crate) fn _is_pipe_constant_heat_flux(&self) -> bool {
        matches!(self, Self::PipeConstantHeatFluxFullyDeveloped)
    }


}

