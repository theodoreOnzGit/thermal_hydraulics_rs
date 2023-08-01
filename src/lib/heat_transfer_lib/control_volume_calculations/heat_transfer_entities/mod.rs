use std::f64::consts::PI;

use uom::si::{f64::*, pressure::atmosphere, power::watt, time::second, length::meter};

use crate::heat_transfer_lib::
thermophysical_properties::{Material, 
    specific_enthalpy::{specific_enthalpy, temperature_from_specific_enthalpy}, density::density, thermal_diffusivity::thermal_diffusivity, specific_heat_capacity::specific_heat_capacity};

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

    /// This vector is meant to house a list of maximum timesteps 
    /// and is meant for auto time stepping 
    pub max_timestep_vector: Vec<Time>,

    /// This vector is meant to house a list of maximum timesteps 
    /// based on conduction only
    pub mesh_stability_lengthscale_vector: Vec<Length>,
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
        let initial_timestep_vector: Vec<Time> = vec![];
        let initial_mesh_stability_lengthscale_vector: Vec<Length> 
        = vec![];

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
            max_timestep_vector:
            initial_timestep_vector,
            mesh_stability_lengthscale_vector:
            initial_mesh_stability_lengthscale_vector,
        }

    }

    #[inline]
    fn advance_timestep(&mut self, timestep: Time) -> Result<(), String>{


        // first thing is to sum up all enthalpy changes
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
        // clear the enthalpy change vector and timestep vector 

        self.rate_enthalpy_change_vector.clear();
        self.max_timestep_vector.clear();
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

        // set time step
        let initial_timestep_vector: Vec<Time> = vec![];
        let mut conduction_stability_lengthscale_vector: 
        Vec<Length> = vec![];

        // if it's a sphere, push the radius to the 
        // conduction_stability_lengthscale_vector

        // we do not discretise along theta or phi 

        conduction_stability_lengthscale_vector.push(ball_radius);



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
                    max_timestep_vector: 
                    initial_timestep_vector,
                    mesh_stability_lengthscale_vector:
                    conduction_stability_lengthscale_vector,
                }
            )
        );


        return Ok(ball_control_vol);
    }

    /// this is a function to determine the relevant time scales 
    /// this one is based on conduction, which calculates timescales 
    /// based on mesh fourier number
    /// for stable conduction
    ///
    /// for an uneven volume, it will just take the shortest of 
    /// these lengthscales to determine a proper time scale
    #[inline]
    fn calculate_conduction_timestep(&self) -> Result<Time,String>{

        // first let us get relevant length scales 
        // for shortest time step, we will get the shortest length 
        // scale

        let lengthscale_stability_vec = 
        self.mesh_stability_lengthscale_vector.clone();


        // initiate a simple loop to find the shortest length scale
        let mut min_lengthscale = Length::new::<meter>(0.0);

        for length in lengthscale_stability_vec.iter() {
            if *length > min_lengthscale {
                min_lengthscale = *length;
            }
        }

        // now we can determine the conduction time scale regardless 
        // of whether it's solid or liquid (no gas implementations yet, 
        // this is alpha code, probably needs a few rewrites)
        //

        let max_mesh_fourier_number: f64 = 0.25;

        let control_vol_material = self.material_control_volume.clone();
        let control_vol_pressure = self.pressure_control_volume.clone();
        let cv_temperature = temperature_from_specific_enthalpy(
            self.material_control_volume, 
            self.current_timestep_control_volume_specific_enthalpy, 
            self.pressure_control_volume)?;


        let thermal_diffusivity_coeff: DiffusionCoefficient = 
        thermal_diffusivity(control_vol_material, 
            cv_temperature, 
            control_vol_pressure)?;

        // now let's obtain the minimum timescale for conduction

        let min_conduction_timescale = max_mesh_fourier_number * 
        min_lengthscale *
        min_lengthscale / 
        thermal_diffusivity_coeff;


        return Ok(min_conduction_timescale);
    }

    /// calculates a time scale based on maximum temperature change 
    /// within one time step, 
    #[inline]
    fn get_max_timestep_based_on_max_temperature_change(
        &self,
        max_temperature_change: TemperatureInterval) 
    -> Result<Time,String> {

        // let's get the net power change first 

        // first thing is to sum up all enthalpy rate of change 
        let mut total_enthalpy_rate_change = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            self.rate_enthalpy_change_vector.clone().iter() {

                total_enthalpy_rate_change += *enthalpy_chg_rate;
            }

        // we have Q = m c_p Delta T
        // Q = power * time 
        //
        // power * time = m c_p  Delta T
        //
        // time = m c_p Delta T / power 
        // for small temperature changes, assume cp constant

        let cv_mass_clone = self.mass_control_volume.clone();
        let cv_material = self.material_control_volume.clone();
        let cv_temperature = self.get_temperature()?;
        let cv_pressure = self.pressure_control_volume.clone();

        let cv_heat_capacity = specific_heat_capacity(
            cv_material,
            cv_temperature,
            cv_pressure,
        )?;

        // now we can calculate time step

        let time_step_based_on_max_temperature_change = 
        cv_mass_clone * cv_heat_capacity * max_temperature_change / 
        total_enthalpy_rate_change;

        // now return it to the environment 

        return Ok(time_step_based_on_max_temperature_change);

    }

    /// compiles a list of time steps based on various criteria, 
    ///
    ///
    #[inline]
    pub fn get_max_timestep(&mut self, 
        max_temperature_change: TemperatureInterval) -> Result<Time, String> {

        // let's calculate conduction time step first, based on 
        // thermal diffusivity and fourier number
        let conduction_timestep = self.calculate_conduction_timestep()?;

        // next I want to push it to the vector 

        self.max_timestep_vector.push(conduction_timestep);

        // let's match the material to solid or liquid 

        let cv_material = self.material_control_volume.clone();

        match cv_material {
            Material::Solid(_) => (),
            Material::Liquid(_) => todo!("need to calculate liquid timestep"),
        }

        // last but not least, check if temperature changes exceed a 
        // certain mount, recommend 0.5 C or 1 C for max temperature 
        // change in one time step
        
        let max_temperature_change_timestep: Time = 
        self.get_max_timestep_based_on_max_temperature_change(
            max_temperature_change)?;

        self.max_timestep_vector.push(max_temperature_change_timestep);


        // initiate a simple loop to find the shortest time scale
        // that is "safe" out of all time steps
        let mut min_timescale = Time::new::<second>(100_f64);

        for time in self.max_timestep_vector.iter() {
            if *time < min_timescale {
                min_timescale = *time;
            }
        }
        return Ok(min_timescale);

    }

}




