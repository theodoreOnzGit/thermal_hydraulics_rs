use uom::si::thermodynamic_temperature::kelvin;
use uom::si::{f64::*, pressure::atmosphere};
use crate::heat_transfer_lib::thermophysical_properties::specific_enthalpy::temperature_from_specific_enthalpy;
use crate::heat_transfer_lib::thermophysical_properties::thermal_conductivity::thermal_conductivity;
use crate::heat_transfer_lib::thermophysical_properties
::{Material, specific_enthalpy::specific_enthalpy};
use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::HeatTransferEntity;

/// for 1D calculations, we need to calculate conductance as well,
/// but there is no area, hence, we have to use a unit area to calculate 
/// the conductance 
pub const UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS: f64 = 1.0;

/// Contains possible heat transfer interactions between the nodes
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum HeatTransferInteractionType {
    /// The user specifies a thermal conductance between the nodes
    /// in units of power/kelvin
    UserSpecifiedThermalConductance(ThermalConductance),

    /// 1D Cartesian Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have two control volumes, each node represents a control 
    /// volume
    ///
    /// // ----------------------------
    /// // |                          |
    /// // *                          *
    /// // |                          |
    /// // ----------------------------
    /// // cv_1                      cv_2
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have one material which determines conductivity 
    /// and then a length which determines the distance between 
    /// the two control volumes
    ///
    SingleCartesianThermalConductanceOneDimension(
        Material,
        XThicknessThermalConduction
    ),

    /// 1D Cartesian Coordinates Thermal Resistance, for solids only
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // cv_1                      cv_2                     cv_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two control volumes
    ///
    /// Information must be passed in as a tuple,
    ///
    ///
    DualCartesianThermalConductance(
        (XThicknessThermalConduction, Length),
        (XThicknessThermalConduction, Length),
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes 
    ///
    /// // -------------------------------------------------------
    /// // |                          |                          |
    /// // *                          *                          *
    /// // |                          |                          |
    /// // -------------------------------------------------------
    /// // cv_1                      cv_2                     cv_3
    ///
    /// between them there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// we have two materials which determines conductivity 
    /// and then two lengths which determines the distance between 
    /// the two control volumes 
    ///
    /// one also needs to determine the 
    /// inner diameter, outer diameter and length of the tube 
    /// 
    ///
    /// 
    ///
    DualCylindricalThermalConductance(
        (Material,RadialCylindricalThicknessThermalConduction),
        (Material,RadialCylindricalThicknessThermalConduction),
        (InnerDiameterThermalConduction, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction
    )
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes along the outer wall
    ///
    /// // ----------------------------
    /// // |                          |                          
    /// // *                          *                          *
    /// // |                          |                         (T_f) 
    /// // ----------------------------
    /// // cv_1                  cv_2                     Fluid_node
    ///
    /// between cv_1 and cv_2 there is a thermal resistance 
    /// based on a q'' = k dT/dx
    ///
    /// between cv_2 and fluid_node, there is convection resistance
    /// specified by a Nusselt Number
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two control volumes
    ///
    ///
    /// For convection, the heat flux from solid surface to fluid 
    /// is:
    ///
    /// q = h A(T_s - T_f)
    /// 
    /// for hA
    /// surface area is calculated by specifying an outer diameter 
    /// and a cylindrical axial length
    ///
    ///
    ///
    CylindricalConductionConvectionThermalConductanceOuterWall(
        (Material,RadialCylindricalThicknessThermalConduction),
        (HeatTransfer, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction),
    ),


    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition(Power),
}

/// For this part, we determine how heat_transfer_entities interact 
/// with each other by using a function 
/// The function will take in a HeatTransferInteractionType enum
/// which you must first initiate
///
/// Then you need to supply two control volumes or more generally 
/// heat_transfer_entities, which can consist of mix of control volumes
/// and boundary conditions
///
/// The function will then calculate the heat transfer between the two 
/// control volumes, and either return a value or mutate the CV objects 
/// using mutable borrows
pub fn link_heat_transfer_entity(entity_1: &mut HeatTransferEntity,
    entity_2: &mut HeatTransferEntity,
    interaction: HeatTransferInteractionType){

    // first thing first, probably want to unpack the enums to obtain 
    // the underlying control volume and BCs

    let control_vol_1: Option<&mut CVType> = match entity_1 {
        HeatTransferEntity::ControlVolume(control_vol_type) 
            => Some(control_vol_type),
        _ => None,
    };

    let control_vol_2: Option<&mut CVType> = match entity_2 {
        HeatTransferEntity::ControlVolume(control_vol_type) 
            => Some(control_vol_type),
        _ => None,
    };

    // I'll pass in both these option types into a function which 
    // calculates specifically for two control volumes


   
}

// the job of this function is to take in a control volume 
// and then mutate it by calculating its interaction
#[inline]
fn calculate_control_volume_serial(
    control_vol_1: &mut CVType,
    control_vol_2: &mut CVType,
    interaction: HeatTransferInteractionType){

    // let me first match my control volumes to their various types
    //
    // // TODO: handle the case of heat addition 
    //
    // // TODO: handle the case of conduction

}

// this function specifically interacts between two single CV nodes
#[inline]
fn caclulate_between_two_singular_cv_nodes(
    single_cv_1: &mut SingleCVNode,
    single_cv_2: &mut SingleCVNode,
    interaction: HeatTransferInteractionType){

    // let's get the two temperatures of the control volumes first
    // so let me get the enthalpies, and then their respective 
    // temperatures 

    let single_cv_1_enthalpy = single_cv_1.
        current_timestep_control_volume_specific_enthalpy.clone();
    let single_cv_2_enthalpy = single_cv_2.
        current_timestep_control_volume_specific_enthalpy.clone();

    // to get the temperatures, we'll need the material as well 
    let single_cv_1_material = single_cv_1.material_control_volume.clone();
    let single_cv_2_material = single_cv_2.material_control_volume.clone();

    // we'll also need to get their pressures 
    let single_cv_1_pressure = single_cv_1.pressure_control_volume.clone();
    let single_cv_2_pressure = single_cv_2.pressure_control_volume.clone();

    // we will now get their respective temperatures 
    let single_cv_1_temperature = temperature_from_specific_enthalpy(
        single_cv_1_material, single_cv_1_enthalpy, single_cv_1_pressure);
    let single_cv_2_temperature = temperature_from_specific_enthalpy(
        single_cv_2_material, single_cv_2_enthalpy, single_cv_2_pressure);

    // now that we got their respective temperatures we can calculate 
    // the thermal conductance between them
    //
    // for conduction for instance, q = kA dT/dx 
    // conductance is watts per kelvin or 
    // q = (kA)/dx * dT
    // conductance here is kA/dx
    // thermal resistance is 1/conductance
    //
    // for convection, we get: 
    // q = hA (Delta T)
    // hA becomes the thermal conductance
    // 
    // // TODO: probably change the unwrap for later
    let thermal_conductance = get_thermal_conductance(
        single_cv_1_temperature.unwrap(), 
        single_cv_2_temperature.unwrap(),
        single_cv_1_pressure, 
        single_cv_2_pressure, 
        interaction);

    // suppose now we have thermal conductance, we can now obtain the 
    // power flow
}

/// this is thermal conductance function 
/// it calculates thermal conductance based on the supplied enum
fn get_thermal_conductance(
    material_temperature_1: ThermodynamicTemperature,
    material_temperature_2: ThermodynamicTemperature,
    material_pressure_1: Pressure,
    material_pressure_2: Pressure,
    interaction: HeatTransferInteractionType) -> Result<ThermalConductance, String> 
{

    let conductance: ThermalConductance = match 
        interaction {
            HeatTransferInteractionType::UserSpecifiedThermalConductance(
                user_specified_conductance) => user_specified_conductance,
            HeatTransferInteractionType
                ::SingleCartesianThermalConductanceOneDimension(
            material,thickness) => get_conductance_single_cartesian_one_dimension(
                    material,
                    material_temperature_1, 
                    material_temperature_2, 
                    material_pressure_1, 
                    material_pressure_2, 
                    thickness)?,
            _ => todo!("define other enum options"),
        };

    return Ok(conductance);
}

/// conductance is watts per kelvin or 
/// q = (kA)/dx * dT
/// conductance here is kA/dx
/// thermal resistance is 1/conductance
fn get_conductance_single_cartesian_one_dimension(
    material: Material,
    material_temperature_1: ThermodynamicTemperature,
    material_temperature_2: ThermodynamicTemperature,
    material_pressure_1: Pressure,
    material_pressure_2: Pressure,
    thickness: XThicknessThermalConduction) -> Result<ThermalConductance,String> 
{

    // note, the question mark here (?) at the end denotes 
    // error propagation.
    //
    // So basically, we do pattern matching. If it's an error,
    // return an error as the result enum 
    // if it's successful, it will obtain the value in the Ok()
    // variant of the enum

    let average_material_pressure: Pressure = 
    (material_pressure_1 + material_pressure_2)/2.0;

    let average_material_temperature_kelvin = 
    material_temperature_1.value + material_temperature_2.value
    -273.15;

    let average_material_temperature = 
    ThermodynamicTemperature::new::<kelvin>(average_material_temperature_kelvin);


    let material_thermal_conductivity = thermal_conductivity(
        material, 
        average_material_temperature, 
        average_material_pressure)?;    

    let cross_sectional_area: Area = Area::new::<uom::si::area::square_meter>(
        UNIT_AREA_SQ_METER_FOR_ONE_DIMENSIONAL_CALCS);

    let material_thickness: Length = thickness.into();

    let material_conductance: ThermalConductance = 
    material_thermal_conductivity * 
    cross_sectional_area / 
    material_thickness;
    

    return Ok(material_conductance);

}


///
/// For this part, we determine how heat_transfer_entities interact 
/// with each other by using a function 
/// The function will take in a HeatTransferInteractionType enum
/// which you must first initiate
///
/// Notes for parallel computation: 
///
/// if we use mutable borrows early, it will be difficult to use 
/// this in a multithreaded manner
///
/// We can use interior mutability to allow for mutex locks and Arc 
/// types to be used so that parallel calculations can be performed 
///
/// lastly, we can just make sure we return a value from the placeholder 
/// function, and then manually map the heat transfer value to each 
/// heat_transfer_entities pair 
fn placeholder_function_suitable_for_parallel_computation(){}
