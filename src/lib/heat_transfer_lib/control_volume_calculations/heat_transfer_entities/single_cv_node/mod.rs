use std::f64::consts::PI;

use uom::si::area::square_meter;
use uom::si::f64::*;

use uom::si::pressure::atmosphere;

use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_interactions::calculations::UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS;
use crate::heat_transfer_lib::
thermophysical_properties::density::density;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::
thermophysical_properties::Material;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;



use super::CVType;
use super::HeatTransferEntity;
use super::InnerDiameterThermalConduction;
use super::OuterDiameterThermalConduction;

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

    /// This vector houses a list of volumetric flowrates coming into 
    /// and out of the control volume 
    /// by convention, positive flowrates mean going into the 
    /// cv, negative flowrates mean flowing out of the cv 
    pub volumetric_flowrate_vector: Vec<VolumeRate>
}

/// here, we have mostly constructors
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
            volumetric_flowrate_vector:
            vec![],
        }

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
                    volumetric_flowrate_vector:
                    vec![],
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



        let one_dimension_cv = 
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
                    volumetric_flowrate_vector:
                    vec![],
                }
            )
        );


        return Ok(one_dimension_cv);
    }

    /// this function constructs cylinder based on length 
    /// dimensions, d and z 
    ///
    /// z is cylinder length
    ///
    /// d (diameter) is 2r (r is radius)
    #[inline]
    pub fn new_cylinder(z: Length, 
        diameter: Length,
        material: Material,
        cv_temperature: ThermodynamicTemperature,
        pressure: Pressure) -> Result<HeatTransferEntity,
    ThermalHydraulicsLibError>{


        // diameter is pi D^2/4
        let cross_sectional_area: Area = diameter * diameter 
        * 0.25 
        * PI;

        let cylinder_vol: Volume = 
        cross_sectional_area * z;

        let cylinder_density: MassDensity = density(
            material,
            cv_temperature,
            pressure)?;

        let one_dimension_mass: Mass = 
        cylinder_density * cylinder_vol;


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

        conduction_stability_lengthscale_vector.push(z);



        let cylinder = 
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
                    volume: cylinder_vol, 
                    max_timestep_vector: 
                    initial_timestep_vector,
                    mesh_stability_lengthscale_vector:
                    conduction_stability_lengthscale_vector,
                    volumetric_flowrate_vector:
                    vec![],
                }
            )
        );


        return Ok(cylinder);
    }

    /// this function constructs a cylindrical shell based on length 
    /// dimensions, id, od and z 
    ///
    /// z is cylinder length
    ///
    /// id (diameter) is 2r (r is radius)
    #[inline]
    pub fn new_cylindrical_shell(z: Length, 
        id: InnerDiameterThermalConduction,
        od: OuterDiameterThermalConduction,
        material: Material,
        cv_temperature: ThermodynamicTemperature,
        pressure: Pressure) -> Result<HeatTransferEntity,
    ThermalHydraulicsLibError>{

        let id: Length = id.into();
        let od: Length = od.into();

        // diameter is pi D^2/4
        let inner_cross_sectional_area: Area = id * id 
        * 0.25 
        * PI;

        let outer_cross_sectional_area: Area = od * od
        * 0.25 
        * PI;
        
        let cross_sectional_area = outer_cross_sectional_area - 
            inner_cross_sectional_area;

        let cylinder_vol: Volume = 
        cross_sectional_area * z;

        let cylinder_density: MassDensity = density(
            material,
            cv_temperature,
            pressure)?;

        let one_dimension_mass: Mass = 
        cylinder_density * cylinder_vol;


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

        conduction_stability_lengthscale_vector.push(z);



        let cylinder = 
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
                    volume: cylinder_vol, 
                    max_timestep_vector: 
                    initial_timestep_vector,
                    mesh_stability_lengthscale_vector:
                    conduction_stability_lengthscale_vector,
                    volumetric_flowrate_vector:
                    vec![],
                }
            )
        );


        return Ok(cylinder);
    }

    /// not all fluid elements are shaped like a cylinder
    ///
    /// fluid may be forced to flow through some odd shaped pipe
    /// This would be defined by a length and a cross sectional area
    #[inline]
    pub fn new_odd_shaped_pipe(z: Length, 
        cross_sectional_area: Area,
        material: Material,
        cv_temperature: ThermodynamicTemperature,
        pressure: Pressure) -> Result<HeatTransferEntity,
    ThermalHydraulicsLibError>{


        let cylinder_vol: Volume = 
        cross_sectional_area * z;

        let cylinder_density: MassDensity = density(
            material,
            cv_temperature,
            pressure)?;

        let one_dimension_mass: Mass = 
        cylinder_density * cylinder_vol;


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

        conduction_stability_lengthscale_vector.push(z);



        let cylinder = 
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
                    volume: cylinder_vol, 
                    max_timestep_vector: 
                    initial_timestep_vector,
                    mesh_stability_lengthscale_vector:
                    conduction_stability_lengthscale_vector,
                    volumetric_flowrate_vector:
                    vec![],
                }
            )
        );


        return Ok(cylinder);
    }
}


/// calculation contains the advance timestep associated function 
///
pub mod calculation;
pub use calculation::*;


/// contains functions to obtain timestep and other things 
pub mod preprocessing;
pub use preprocessing::*;

