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

    /// pipe nusselt number using custom Gnielinski correlation 
    /// for laminar, turbulent and transition region 
    ///
    /// laminar flow assumes constant heat flux  (Nu = 4.354)
    ///
    /// Correlation be like:
    /// Nu = C (Re^m - 280.0) Pr_f^0.4 ( 1.0 + (D_e/l)^(2/3) ) ( Pr_f / Pr_w )^0.25
    /// User must supply C and m 
    ///
    /// For low Re flows, Nu = 4.36 is used. 
    /// Turbulence is assumed at Re = 4000 taking pipes as a reference.
    /// The transition regime is around Re = 2300 - 4000 
    /// this is taken from the Re for transition in pipes 
    /// THESE MAY NOT BE APPLICABLE IN THIS CASE
    CustomGnielinskiGeneric(GnielinskiData, Ratio, f64),

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

    /// Ideal 1e9 
    /// Just returns a Nusselt number of 10^9 
    /// which may be suitable as an approximation for heat exchangers 
    IdealNusseltOneBillion,
}

impl NusseltCorrelation {


    /// gets the nusselt based on user choice of of correlation
    pub fn try_get(&self) -> Result<Ratio, ThermalHydraulicsLibError> {
        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric(data) => {
                return data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::CustomGnielinskiGeneric(
                data, correlation_coefficient_c, reynolds_exponent_m
            ) => 
            {
                return data.get_nusselt_for_custom_developing_flow
                    (*correlation_coefficient_c,*reynolds_exponent_m)
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
            NusseltCorrelation::IdealNusseltOneBillion => {
                Ratio::new::<ratio>(1e9_f64)
            },
        };

        return Ok(nusselt_number);
    }

    /// gets an estimate for the nusselt number based on user choice 
    /// of correlation, ignores wall temperature 
    ///
    /// note that this uses clone, so it's quite resource heavy
    pub fn estimate_based_on_prandtl_and_reynolds_no_wall_correction(&self,
    bulk_prandtl_number_input: Ratio,
    reynolds_number_input: Ratio,) -> Result<Ratio, ThermalHydraulicsLibError>{

        let nusselt_number: Ratio = 
        match self {
            NusseltCorrelation::PipeGnielinskiGeneric(data) => {

                let mut modified_data = data.clone();
                modified_data.prandtl_wall = bulk_prandtl_number_input;
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_developing_flow();
            },
            NusseltCorrelation::CustomGnielinskiGeneric(
                data, correlation_coefficient_c, reynolds_exponent_m
            ) => 
            {
                let mut modified_data = data.clone();
                modified_data.prandtl_wall = bulk_prandtl_number_input;
                modified_data.prandtl_bulk = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_custom_developing_flow
                    (*correlation_coefficient_c,*reynolds_exponent_m)
            },
            NusseltCorrelation::PipeGnielinskiTurbulent(data) => {
                let mut modified_data = data.clone();
                modified_data.prandtl_wall = bulk_prandtl_number_input;
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
                modified_data.prandtl_wall = bulk_prandtl_number_input;
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
                modified_data.prandtl_wall = bulk_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;

                return modified_data.ciet_version_2_heater_prandtl_corrected();
            },
            NusseltCorrelation::IdealNusseltOneBillion => {
                Ratio::new::<ratio>(1e9_f64)
            },
        };

        return Ok(nusselt_number);
    }

    /// gets an estimate for the nusselt number based on user choice 
    /// of correlation, includes wall temperature correction
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
            NusseltCorrelation::CustomGnielinskiGeneric(
                data, correlation_coefficient_c, reynolds_exponent_m
            ) => 
            {
                let mut modified_data = data.clone();
                modified_data.prandtl_wall = bulk_prandtl_number_input;
                modified_data.prandtl_bulk = wall_prandtl_number_input;
                modified_data.reynolds = reynolds_number_input;
                return modified_data.get_nusselt_for_custom_developing_flow
                    (*correlation_coefficient_c,*reynolds_exponent_m)
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
            NusseltCorrelation::IdealNusseltOneBillion => {
                Ratio::new::<ratio>(1e9_f64)
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

