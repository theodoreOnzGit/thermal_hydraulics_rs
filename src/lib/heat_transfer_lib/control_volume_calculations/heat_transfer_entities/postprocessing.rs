use super::{HeatTransferEntity, CVType, BCType, ArrayCVType};

use uom::si::{f64::*, pressure::atmosphere, power::watt, time::second, length::meter};

impl HeatTransferEntity {

    /// gets the temperature of the HeatTransferEntity 
    /// usually control volume at the current timestep
    pub fn temperature(entity: &mut HeatTransferEntity) -> 
    Result<ThermodynamicTemperature, String> {

        let cv_temperature_result = match entity {
            Self::ControlVolume(cv_type) => {
                match cv_type {
                    CVType::SingleCV(single_cv) => {
                        single_cv.get_temperature()
                    },
                    CVType::ArrayCV(cv) => {
                        cv.get_bulk_temperature()
                    },
                }
            },
            Self::BoundaryConditions(bc_type) => {
                match bc_type {
                    BCType::UserSpecifiedTemperature(temperature) => {
                        Ok(*temperature)
                    },
                    BCType::UserSpecifiedHeatFlux(_) | 
                        BCType::UserSpecifiedHeatAddition(_) => {
                            return Err("getting temperature not \n 
                                implemented for BoundaryConditions".to_string())
                        },
                }
            },
        };


        return cv_temperature_result;
    }

    /// gets a vector of temperatures
    pub fn temperature_vector(entity: &mut HeatTransferEntity) ->
    Result<Vec<ThermodynamicTemperature>,String> {


        // oof nested matching, kind of ugly but I'll live with it for now
        let cv_temperature_vector_result: 
        Result<Vec<ThermodynamicTemperature>, String>
        = match entity {
            HeatTransferEntity::ControlVolume(cv_type) => {
                let temp_vector_result:
                Result<Vec<ThermodynamicTemperature>, String>
                = match cv_type{
                    CVType::SingleCV(single_cv) => {
                        // the single cv only has one temperature anyway
                        let temperature = single_cv.get_temperature()?;

                        let mut temp_vector: Vec<ThermodynamicTemperature>
                        = vec![];

                        temp_vector.push(temperature);

                        Ok(temp_vector)
                    },
                    CVType::ArrayCV(array_cv_type) => {
                        let temp_vector_result:
                        Result<Vec<ThermodynamicTemperature>, String>;

                        temp_vector_result = match array_cv_type {
                            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                                cartesian_array_cv.get_temperature_vector()
                            },
                        };
                        temp_vector_result
                    },
                };
                temp_vector_result
            },
            Self::BoundaryConditions(bc_type) => {
                match bc_type {
                    BCType::UserSpecifiedTemperature(temperature) => {

                        let mut temp_vector: Vec<ThermodynamicTemperature>
                        = vec![];

                        temp_vector.push(*temperature);
                        Ok(temp_vector)
                    },
                    BCType::UserSpecifiedHeatFlux(_) | 
                        BCType::UserSpecifiedHeatAddition(_) => {
                            return Err("getting temperature not \n 
                                implemented for BoundaryConditions".to_string())
                        },
                }
            },
        };

        cv_temperature_vector_result
    }
}
