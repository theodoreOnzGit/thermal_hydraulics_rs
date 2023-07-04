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


/// Contains possible heat transfer interactions between the nodes
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum HeatTransferInteractionTypes {
    /// The user specifies a thermal conductance between the nodes
    /// in units of power/kelvin
    UserSpecifiedThermalConductance(ThermalConductance),

    /// 1D Cartesian Coordinates Thermal Resistance
    ///
    /// basically have two nodes 
    ///
    /// // ----------------------------
    /// // |                          |
    /// // *                          *
    /// // |                          |
    /// // ----------------------------
    /// // node_1                  node_2
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have one material which determines conductivity 
    /// and then a length which determines the distance between 
    /// the two nodes
    ///
    SingleCartesianThermalResistance(Material,XThicknessThermalConduction),

    /// 1D Cartesian Coordinates Thermal Resistance
    ///
    /// basically have three nodes 
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // node_1                  node_2                     node_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two nodes
    ///
    /// Information must be passed in as a tuple,
    ///
    ///
    DualCartesianThermalResistance(
        (XThicknessThermalConduction, Length),
        (XThicknessThermalConduction, Length),
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    ///
    /// basically have three nodes 
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // node_1                  node_2                     node_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two nodes 
    ///
    /// one also needs to determine the 
    /// inner diameter, outer diameter and length of the tube 
    /// 
    /// // TODO: 
    /// This is the bare minimum we need, 
    /// Though TBH, I also want the compiler to alert the user 
    /// as to what kind of stuff to put in
    ///
    /// Like a type wrapper,
    /// I'll do this after lunch (TBD)
    ///
    DualCylindricalThermalResistance(
        (Material,Length),
        (Material,Length),
        (Length, Length, Length)
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    ///
    /// basically have three nodes along the outer wall
    ///
    /// // ----------------------------
    /// // |                          |                          
    /// // *                          *                          *
    /// // |                          |                         (T_f) 
    /// // ----------------------------
    /// // node_1                  node_2                     Fluid_node
    ///
    /// between node_1 and node_2 there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// between node_2 and fluid_node, there is convection resistance
    /// specified by a Nusselt Number
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two nodes
    ///
    ///
    /// For convection, the heat flux from solid surface to fluid 
    /// is:
    ///
    /// q'' = h (T_s - T_f)
    ///
    ///
    CylindricalConductionConvectionThermalResistanceOuterWall(
        Material,
        Length,
        Material,
        Length
    ),


    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}


/// XThicknessThermalConduction is essentially a struct containing 
/// one length describing a thickness in cartesian coordinates 
/// for thermal conduction.
///
/// This type is meant for an input for various enums and functions 
/// in this crate.
///
/// It is meant to guide the user so that they know what the 
/// length inputs represents.
///
/// ```rust
///
/// use uom::si::length::meter;
/// use uom::si::f64::*;
/// use thermal_hydraulics_rs::heat_transfer_lib::
/// control_volume_calculations:: heat_transfer_entities:: 
/// XThicknessThermalConduction;
///
/// // let's say you have a thickness of 0.5 which you want to describe
///
/// let thickness_of_wall = Length::new::<meter>(0.5);
///
/// // we first need to convert it into an XThicknessThermalConduction
/// // type first
/// let wall_thickness_input = XThicknessThermalConduction::from(
/// thickness_of_wall);
///
/// // to convert it back into a length type, we use the into() method 
///
/// let thickness_wall_for_calculation: Length = 
/// wall_thickness_input.into();
/// 
/// // both these are the same
/// assert_eq!(thickness_of_wall, thickness_wall_for_calculation);
/// ```
/// 
/// 
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct XThicknessThermalConduction {
    thickness: Length,
}

impl From<Length> for XThicknessThermalConduction {
    fn from(thickness: Length) -> Self{
        Self { thickness }
    }
}

impl Into<Length> for XThicknessThermalConduction {
    fn into(self) -> Length {
        self.thickness
    }
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


