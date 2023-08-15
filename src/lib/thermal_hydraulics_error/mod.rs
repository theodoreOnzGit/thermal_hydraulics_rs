use thiserror::Error;

/// Master Error type of this crate
#[derive(Debug, Error)]
pub enum ThermalHydraulicsError {
    /// linear algebra error
    #[error("linear algebra error")]
    LinalgError(#[from] ndarray_linalg::error::LinalgError),

    /// empty mass flowrate vector error 
    ///
    /// this case is where the mass flowrate vector in a control 
    /// volume is empty, 
    /// so we can't calculate a courant number

    #[error("cannot calculate courant number: mass flowrate \n 
        there is no mass flows going in or out of your \n 
        control volume")]
    CourantMassFlowVectorEmpty,


    /// it's a generic error which
    #[error("Placeholder Error Type for Strings{0} ")]
    GenericStringError(String),
}

/// 
impl From<String> for ThermalHydraulicsError {
    fn from(value: String) -> Self {
        Self::GenericStringError(value)
    }
}

impl Into<String> for ThermalHydraulicsError {
    fn into(self) -> String {
        match self {
            ThermalHydraulicsError::LinalgError(_) => {
                self.to_string()
            },
            ThermalHydraulicsError::CourantMassFlowVectorEmpty => {
                self.to_string()
            },
            ThermalHydraulicsError::GenericStringError(string) => {
                string
            },
        }
    }
}

