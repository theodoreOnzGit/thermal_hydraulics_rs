use std::f64::consts::PI;

use uom::si::{f64::*, pressure::atmosphere, power::watt};

use crate::heat_transfer_lib::
thermophysical_properties::{Material, 
    specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy}, density::density};

/// Contains entities which transfer heat and interact with each 
/// other
///
/// for example, control volumes and boundary conditions
#[derive(Debug,Clone,PartialEq)]
pub enum HeatTransferEntity {
    /// Contains a list of ControlVolumeTypes
    ControlVolume(CVType),
    /// Contains a list of Boundary conditions
    BoundaryConditions(BCType)
}
/// placeholder, I'd like to have some associated functions to 
/// deal with the HeatTransferEntity type
///
/// probably one to get the courant number, 
/// and second, to use a timestep to calculate the new enthalpy 
/// and update enthalpy
/// 
/// last but not least, extract temperatures for sensing purposes
impl HeatTransferEntity {

    /// for control volumes, this method allows you to 
    /// calculate the enthalpy of the next timestep and 
    /// set the 
    /// current timestep enthalpy as the enthalpy calculated  
    /// for the next timestep
    ///
    /// you are required to explicitly provide a timestep for this 
    pub fn advance_timestep(entity: &mut HeatTransferEntity,
    timestep: Time) -> Result<(), String> {

        // first match CV or BC, 
        // Boundary conditions don't need to advance timestep
        // so we can leave them be (it should return an Ok(()) value 
        // rather than an Err() value)

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => return Ok(()),
        };

        // once I have the cv_type enum, match it again

        let cv_advance_result: Result<(), String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.advance_timestep(timestep)
                },
                CVType::ArrayCV => return Err("not implemented".to_string()),
            };

        return cv_advance_result;
    }

    /// gets the temperature of the HeatTransferEntity 
    /// usually control volume at the current timestep
    pub fn temperature(entity: &mut HeatTransferEntity) -> 
    Result<ThermodynamicTemperature, String> {

        let control_vol_type = match entity {
            Self::ControlVolume(cv_type) => cv_type,
            Self::BoundaryConditions(_) => 
                return Err("getting temperature not \n 
                    implemented for BoundaryConditions".to_string()),
        };

        // once I have the cv_type enum, match it again

        let cv_temperature_result: 
        Result<ThermodynamicTemperature, String> = match 
            control_vol_type {
                CVType::SingleCV(single_cv) => {
                    single_cv.get_temperature()
                },
                CVType::ArrayCV => return Err("not implemented".to_string()),
            };

        return cv_temperature_result;
    }

}

/// To determine heat transfer between two control volumes or 
/// generally, two heat transfer entities, one must determine 
/// which control volume is in front and which is at the back
/// 
/// This type would tell the solver that this control volume is 
/// in the front
/// 
#[derive(Debug,Clone,PartialEq)]
pub struct FrontHeatTransferEntity {
    entity: HeatTransferEntity,
}

impl From<HeatTransferEntity> for FrontHeatTransferEntity {
    fn from(entity: HeatTransferEntity) -> Self{
        Self { entity }
    }
}

impl Into<HeatTransferEntity> for FrontHeatTransferEntity {
    fn into(self) -> HeatTransferEntity {
        self.entity
    }
}

/// To determine heat transfer between two control volumes or 
/// generally, two heat transfer entities, one must determine 
/// which control volume is in front and which is at the back
/// 
/// This type would tell the solver that this control volume is 
/// in the back
/// 
#[derive(Debug,Clone,PartialEq)]
pub struct BackHeatTransferEntity {
    entity: HeatTransferEntity,
}

impl From<HeatTransferEntity> for BackHeatTransferEntity {
    fn from(entity: HeatTransferEntity) -> Self{
        Self { entity }
    }
}

impl Into<HeatTransferEntity> for BackHeatTransferEntity {
    fn into(self) -> HeatTransferEntity {
        self.entity
    }
}

/// Contains Types of Control Volumes (CVs)
#[derive(Debug,Clone,PartialEq)]
pub enum CVType {
    /// This CV is the most basic,  it can be represented by a single 
    /// point or node
    SingleCV(SingleCVNode),
    /// Array CVs are collections of SingleCVs, but do not require the 
    /// user to manually specify the connections between the SingleCVs
    ArrayCV,
}

/// Contains all the types of Boundary Conditions (BCs) you can use 
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum BCType {
    /// The user specifies a fixed temperature for the BC
    UserSpecifiedTemperature(ThermodynamicTemperature),
    /// The user specifies a heat flux for the BC
    /// the uom type is heat flux density in power/area
    UserSpecifiedHeatFlux(HeatFluxDensity),
    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}



/// Contains different length types for use in defining interactions 
/// between heat transfer entities
///
/// This is to make it clear to the user exactly what 
/// kind of length we are specifying 
///
/// For example, to define an annular hollow cylinder, we need three 
/// lengths:
///
/// the zLength of the cylinder,
/// inner diameter 
/// outer diameter
///
/// each of these have type Length. The user would have to read 
/// the documentation as to what kind of length is being 
/// specified by the user 
///
/// Hence, I'm making some extra types to force the compiler to tell you 
/// (the user) what kind of length you need to specify
///
/// 
///
pub mod heat_transfer_dimensions;
pub use heat_transfer_dimensions::*;


/// SingleCVNode (single control volume node) represents 
/// the control volume with a fixed point
///
/// The idea for a SingleCVNode, is for it to contain information 
/// about a control volume. 
///
/// One can then connect these control volumes with other control 
/// volumes and then specify the interaction or heat transfer between
/// adjacent Control Volumes CVs and Boundary Conditions BCs
///
/// The Control Volume is initiated with a temperature and material 
/// type, this would help determine the control volume's specific 
/// energy,
/// the mass of the system must also be specified
///
/// The changes can be pushed to a vector called the enthalpy 
/// change vector
///
/// At the end of the timestep, the next_timestep_specific_enthalpy 
/// is calculated by the current_timestep_control_volume_specific_enthalpy
/// plus the enthalpy changes in the vector
///
/// The temperature can then be calculated from the 
/// next_timestep_specific_enthalpy
///
/// 
///
#[derive(Debug,Clone,PartialEq)]
pub struct SingleCVNode {

    /// specific enthalpy at present timestep, set using 
    /// the temperature and material type
    pub current_timestep_control_volume_specific_enthalpy: AvailableEnergy,
    /// specific enthalpy at next timestep, used to calculate 
    /// temperature
    pub next_timestep_specific_enthalpy: AvailableEnergy,

    /// contains rate of change of the specific enthalpy due to changes
    /// once courant_number is determined, we would use the correct 
    /// timestep to multiply the power into an overall enthalpy change
    pub rate_enthalpy_change_vector: Vec<Power>,

    /// control volume mass 
    pub mass_control_volume: Mass,

    /// control volume material 
    pub material_control_volume: Material,

    /// control volume pressure 
    pub pressure_control_volume: Pressure,

    /// volume of the control volume 
    pub volume: Volume,
}

impl SingleCVNode {
    /// to initiate the control volume, use this constructor,
    /// which means we supply the temperature, material type 
    /// and mass of the CV
    ///
    /// assumes properties are at atmospheric pressure
    pub fn new(cv_temperature: ThermodynamicTemperature,
        cv_material: Material,
        cv_mass: Mass,
        cv_volume: Volume) -> SingleCVNode {

        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let cv_enthalpy: AvailableEnergy = 
        match specific_enthalpy(
            cv_material, 
            cv_temperature, 
            atmospheric_pressure) {
                Ok(specific_enthalpy) => specific_enthalpy,
                Err(error_msg) => panic!("{}", error_msg),

        };

        // we store a vector of power at current timestep,
        // once courant_number is decided, then we shall decide the 
        // timestep
        let cv_rate_enthalpy_change_vec: Vec<Power> = vec![];

        return Self{
            current_timestep_control_volume_specific_enthalpy : 
            cv_enthalpy,
            next_timestep_specific_enthalpy : 
            cv_enthalpy,
            rate_enthalpy_change_vector : 
            cv_rate_enthalpy_change_vec,
            mass_control_volume : 
            cv_mass,
            material_control_volume: 
            cv_material,
            pressure_control_volume:
            atmospheric_pressure,
            volume: cv_volume,
        }

    }

    #[inline]
    fn advance_timestep(&mut self, timestep: Time) -> Result<(), String>{


        let mut total_enthalpy_rate_change = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            self.rate_enthalpy_change_vector.clone().iter() {

                total_enthalpy_rate_change += *enthalpy_chg_rate;
            }

        let enthalpy_next_timestep = total_enthalpy_rate_change * 
        timestep.clone() +
        self.current_timestep_control_volume_specific_enthalpy.
            clone()* self.mass_control_volume.clone();

        let specific_enthalpy_next_timestep = 
        enthalpy_next_timestep/self.mass_control_volume.clone();


        self.next_timestep_specific_enthalpy 
            = specific_enthalpy_next_timestep;

        // at the end of each timestep, set 
        // current_timestep_control_volume_specific_enthalpy
        // to that of the next timestep

        self.current_timestep_control_volume_specific_enthalpy
            = specific_enthalpy_next_timestep;

        // now things get a bit clunky here, I will need to set 
        // the mass of the control volume by the temperature after 
        // I calculated it using control volume mass from the last 
        // timestep. For large changes in density, this may not 
        // work well (instabilities??) but for liquids, I hope it's 
        // okay

        self.set_liquid_cv_mass_from_temperature()?;
        // clear the enthalpy change vector

        self.rate_enthalpy_change_vector.clear();
        // increase timestep (last step)

        return Ok(());
    }

    #[inline]
    fn get_temperature(&self) -> 
    Result<ThermodynamicTemperature, String>{

        let cv_temperature = temperature_from_specific_enthalpy(
            self.material_control_volume, 
            self.current_timestep_control_volume_specific_enthalpy, 
            self.pressure_control_volume);

        return cv_temperature;
    }

    /// this function takes the temperature of the control volume 
    /// to find its density and set its mass
    ///
    /// only changes mass for liquids, not solids
    #[inline]
    fn set_liquid_cv_mass_from_temperature(&mut self) 
    -> Result<(), String>{

        let cv_temperature = temperature_from_specific_enthalpy(
            self.material_control_volume, 
            self.current_timestep_control_volume_specific_enthalpy, 
            self.pressure_control_volume)?;

        let cv_density = density(
            self.material_control_volume, 
            cv_temperature, 
            self.pressure_control_volume)?;

        // we are going to need a cv_volume
        let cv_volume = self.volume;

        let new_cv_mass = cv_density * cv_volume;

        match self.material_control_volume {
            // for solids, do not set anything
            Material::Solid(_) => return Ok(()),
            Material::Liquid(_) => {
                self.mass_control_volume = new_cv_mass;
                return Ok(());
            },
        }

    }

    /// this function constructs control volume based on spherical 
    /// dimensions,
    /// it will return a HeatTransferEntity so that you don't have 
    /// to do all the packaging manually
    #[inline]
    pub fn new_sphere(diameter: Length, 
        material: Material,
        cv_temperature: ThermodynamicTemperature,
        pressure: Pressure) -> Result<HeatTransferEntity,String>{


        let ball_radius: Length = diameter * 0.5;

        let ball_volume: Volume = 4.0/3.0 * PI * 
        ball_radius * ball_radius * ball_radius;

        let ball_density: MassDensity = density(
            material,
            cv_temperature,
            pressure)?;

        let ball_mass: Mass = ball_density * ball_volume;


        let enthalpy = specific_enthalpy(
            material, 
            cv_temperature, 
            pressure)?;

        let ball_control_vol = 
        HeatTransferEntity::ControlVolume(
            CVType::SingleCV(
                    SingleCVNode { 
                    current_timestep_control_volume_specific_enthalpy: 
                    enthalpy, 
                    next_timestep_specific_enthalpy: 
                    enthalpy, 
                    rate_enthalpy_change_vector: 
                    vec![], 
                    mass_control_volume: ball_mass, 
                    material_control_volume: material, 
                    pressure_control_volume: pressure,
                    volume: ball_volume, 
                }
            )
        );


        return Ok(ball_control_vol);
    }

}




