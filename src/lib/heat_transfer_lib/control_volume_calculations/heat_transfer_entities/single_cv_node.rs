use std::f64::consts::PI;

use uom::si::area::square_meter;
use uom::si::f64::*;
use uom::si::length::meter;
use uom::si::power::watt;
use uom::si::pressure::atmosphere;
use uom::si::time::second;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::calculations::UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS;
use crate::heat_transfer_lib::
thermophysical_properties::specific_heat_capacity::specific_heat_capacity;
use crate::heat_transfer_lib::
thermophysical_properties::thermal_diffusivity::thermal_diffusivity;
use crate::heat_transfer_lib::
thermophysical_properties::density::density;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::Material;



use super::CVType;
use super::HeatTransferEntity;

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
#[derive(Debug,Clone,PartialEq,Default)]
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

    /// this function performs necessary calculations to move 
    /// the state of the control volume to the next time step
    ///
    /// calculates the new enthalpy of the 
    /// and cleans out the all power vectors and time step vectors
    #[inline]
    pub fn advance_timestep(&mut self, timestep: Time) -> Result<(), String>{


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

    /// gets the temperature of the control volume at the 
    /// CURRENT timestep
    #[inline]
    pub fn get_temperature(&self) -> 
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
    pub (in crate::heat_transfer_lib::control_volume_calculations)
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

    /// this function constructs 1d control volume based on spherical 
    /// dimensions,
    /// it will return a HeatTransferEntity so that you don't have 
    /// to do all the packaging manually
    #[inline]
    pub fn new_one_dimension_volume(length: Length, 
        material: Material,
        cv_temperature: ThermodynamicTemperature,
        pressure: Pressure) -> Result<HeatTransferEntity,String>{


        let basis_area: Area = Area::new::<square_meter>( 
            UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS);

        let one_dimension_volume: Volume = 
        basis_area * length;

        let one_dimension_density: MassDensity = density(
            material,
            cv_temperature,
            pressure)?;

        let one_dimension_mass: Mass = 
        one_dimension_density * one_dimension_volume;


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

        conduction_stability_lengthscale_vector.push(length);



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
                    mass_control_volume: one_dimension_mass, 
                    material_control_volume: material, 
                    pressure_control_volume: pressure,
                    volume: one_dimension_volume, 
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
    pub fn calculate_conduction_timestep(&self) -> Result<Time,String>{

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

    /// appends timestep constrained to fourier number stability
    #[inline]
    pub fn append_conduction_mesh_stability_timestep(&mut self, 
        lengthscale: Length) -> Result<Time, String> {

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
        lengthscale *
        lengthscale / 
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
        //
        // inspired by GeN-Foam and other papers.
        
        let max_temperature_change_timestep: Time = 
        self.get_max_timestep_based_on_max_temperature_change(
            max_temperature_change)?;

        self.max_timestep_vector.push(max_temperature_change_timestep);

        // we still need the timestep changes from the control volumes 
        // linked to this cv.


        // initiate a simple loop to find the shortest time scale
        // that is "safe" out of all time steps
        let mut min_timescale = Time::new::<second>(75_f64);

        for time in self.max_timestep_vector.iter() {
            if *time < min_timescale {
                min_timescale = *time;
            }
        }
        return Ok(min_timescale);

    }

}




