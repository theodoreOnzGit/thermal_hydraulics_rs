use crate::{thermal_hydraulics_error::ThermalHydraulicsLibError, heat_transfer_lib::thermophysical_properties::density::try_get_rho};

use super::HeatTransferEntity;
use super::BCType;
use super::CVType;


use uom::si::f64::*;

impl HeatTransferEntity {

    #[inline]
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
                        Ok(cv.get_bulk_temperature().unwrap())
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

    /// gets bulk temperature of the heat transfer entity
    #[inline]
    pub fn get_bulk_temperature(&mut self) -> Result<ThermodynamicTemperature,
    ThermalHydraulicsLibError> {

        let bulk_temp: ThermodynamicTemperature = 
        HeatTransferEntity::temperature(self)?;

        Ok(bulk_temp)
    }

    /// gets a vector of temperatures
    #[inline]
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


                        let temp_vector_result = 
                        array_cv_type.get_temperature_vector();


                        match temp_vector_result {
                            Ok(vec) => Ok(vec),
                            Err(err) => Err(err.to_string()),
                        }
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

    /// gets temperature vector of this HeatTransferEntity 
    #[inline]
    pub fn get_temperature_vector(&mut self) ->
    Result<ThermodynamicTemperature, ThermalHydraulicsLibError>{

        let temp_vec: ThermodynamicTemperature = 
        HeatTransferEntity::temperature(self)?;

        Ok(temp_vec)
    }

    /// density vector 
    /// attempts to get a vector of densities
    #[inline]
    pub fn density_vector(entity: &mut HeatTransferEntity) ->
    Result<Vec<MassDensity>,ThermalHydraulicsLibError> {

        let temperature_vector = Self::temperature_vector(entity)?;

        let material = match entity {
            HeatTransferEntity::ControlVolume(cv) => {
                match cv {
                    CVType::SingleCV(single_cv) => {
                        single_cv.material_control_volume.clone()
                    },
                    CVType::ArrayCV(array_cv) => {
                        array_cv.get_array_cv_material().unwrap()
                    },
                }
            },
            HeatTransferEntity::BoundaryConditions(_) => {
                return Err(ThermalHydraulicsLibError::
                    NotImplementedForBoundaryConditions(
                        "density not implemented for BC".to_string()
                    ));
            },
        };

        let pressure = match entity {
            HeatTransferEntity::ControlVolume(cv) => {
                match cv {
                    CVType::SingleCV(single_cv) => {
                        single_cv.pressure_control_volume
                    },
                    CVType::ArrayCV(array_cv) => {
                        array_cv.get_array_cv_pressure().unwrap()
                    },
                }
            },
            HeatTransferEntity::BoundaryConditions(_) => {
                return Err(ThermalHydraulicsLibError::
                    NotImplementedForBoundaryConditions(
                        "density not implemented for BC".to_string()
                    ));
            },

        };

        // get density 

        let mut density_vector: Vec<MassDensity> = vec![];

        for temperature in temperature_vector {
            let density = try_get_rho(
                material,
                temperature,
                pressure
            )?;

            density_vector.push(density);

        }

        // return density vector
        Ok(density_vector)


    }
}
