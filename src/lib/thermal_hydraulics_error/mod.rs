use thiserror::Error;
use tuas_boussinesq_solver::tuas_lib_error::TuasLibError;

/// Master Error type of this crate
#[derive(Debug, Error)]
pub enum ThermalHydraulicsLibError {
    /// linear algebra error
    #[error("linear algebra error")]
    LinalgError(#[from] ndarray_linalg::error::LinalgError),

    /// it's a generic error which is a placeholder since I used 
    /// so many string errors
    #[error("Placeholder Error Type for Strings{0} ")]
    GenericStringError(String),


    /// TUAS Boussinesq Solver error 
    #[error("Tuas Boussinesq Solver Error")]
    TUASError(TuasLibError),
    
}

///  converts ThermalHydraulicsLibError from string error
impl From<String> for ThermalHydraulicsLibError {
    fn from(value: String) -> Self {
        Self::GenericStringError(value)
    }
}

impl Into<String> for ThermalHydraulicsLibError {
    fn into(self) -> String {
        match self {
            ThermalHydraulicsLibError::LinalgError(_) => {
                self.to_string()
            },
            ThermalHydraulicsLibError::GenericStringError(string) => {
                string
            },
            ThermalHydraulicsLibError::TUASError(_) => {
                self.to_string()
            },

        }
    }
}

