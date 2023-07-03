use uom::si::{f64::*, pressure::atmosphere};

use crate::heat_transfer_lib::thermophysical_properties::{Material, specific_enthalpy::specific_enthalpy};

/// Contains entities which transfer heat and interact with each 
/// other
///
/// for example, control volumes and boundary conditions
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum HeatTransferEntities {
    /// Contains a list of ControlVolumeTypes
    ControlVolume(CVTypes),
    /// Contains a list of selectable liquids
    BoundaryConditions(BCTypes)
}


/// Contains Types of Control Volumes (CVs)
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum CVTypes {
    /// This CV is the most basic,  it can be represented by a single 
    /// point or node
    SingleCV,
    /// Array CVs are collections of SingleCVs, but do not require the 
    /// user to manually specify the connections between the SingleCVs
    ArrayCV,
}

/// Contains all the types of Boundary Conditions (BCs) you can use 
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum BCTypes {
    /// The user specifies a fixed temperature for the BC
    UserSpecifiedTemperature(ThermodynamicTemperature),
    /// The user specifies a heat flux for the BC
    /// the uom type is heat flux density in power/area
    UserSpecifiedHeatFlux(HeatFluxDensity),
    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}


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
    current_timestep_control_volume_specific_enthalpy: AvailableEnergy,
    /// specific enthalpy at next timestep, used to calculate 
    /// temperature
    next_timestep_specific_enthalpy: AvailableEnergy,

    /// contains changes to the specific enthalpy due to changes
    enthalpy_change_vector: Vec<Energy>,

    /// control volume mass 
    mass_control_volume: Mass,
}

impl SingleCVNode {
    /// to initiate the control volume, use this constructor,
    /// which means we supply the temperature, material type 
    /// and mass of the CV
    ///
    /// assumes proeprties are at atmospheric pressure
    pub fn new(cv_temperature: ThermodynamicTemperature,
        cv_material: Material,
        cv_mass: Mass) -> SingleCVNode {

        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

        let cv_enthalpy: AvailableEnergy = 
        match specific_enthalpy(
            cv_material, 
            cv_temperature, 
            atmospheric_pressure) {
                Ok(specific_enthalpy) => specific_enthalpy,
                Err(error_msg) => panic!("{}", error_msg),

        };

        let cv_enthalpy_change_vec: Vec<Energy> = vec![];

        return Self{
            current_timestep_control_volume_specific_enthalpy : 
            cv_enthalpy,
            next_timestep_specific_enthalpy : 
            cv_enthalpy,
            enthalpy_change_vector : 
            cv_enthalpy_change_vec,
            mass_control_volume : 
            cv_mass,
        }

    }
}


