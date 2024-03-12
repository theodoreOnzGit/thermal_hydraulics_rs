use crate::boussinesq_solver::control_volume_dimensions::*;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::Material;
use uom::si::f64::*;

use super::heat_transfer_geometry::DataDualCartesianThermalConductanceThreeDimension;


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


    /// suppose there are two blocks with the same cross sectional 
    /// area, each of its own thickness and material makeup 
    ///
    /// this is DualCartesianThermalConductanceThreeDimension
    /// we have three dimensional blocks, but the conduction is along 
    /// the thickness of the block, tube or cylinder
    DualCartesianThermalConductanceThreeDimension(
        DataDualCartesianThermalConductanceThreeDimension
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
        (Material, XThicknessThermalConduction),
        (Material, XThicknessThermalConduction),
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
    /// the first material and thickness argument represents 
    /// cv_1 to cv_2 (the inner shell)
    ///
    /// and the second entry pertains to the outer shell 
    /// cv_2 to cv_3, or the outer shell
    ///
    /// 
    ///
    DualCylindricalThermalConductance(
        (Material,RadialCylindricalThicknessThermalConduction),
        (Material,RadialCylindricalThicknessThermalConduction),
        (InnerDiameterThermalConduction, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction)
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes along the outer wall
    ///
    /// -------------------------------------------------------> r
    /// // ----------------------------
    /// // |                          |                          
    /// // * solid_cv_1               *                          *
    /// // |                          |                         (T_f) 
    /// // ----------------------------
    /// //                        solid_surface              Fluid_node
    ///
    /// Where r is the radius 
    /// basically the liquid is on the outside (larger r)
    ///
    /// between solid_cv_1 and the solid_surface 
    /// cv_2 there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// between solid_surface and fluid_node, there is convection resistance
    /// specified by a Nusselt Number so that we get a heat transfer 
    /// coefficient
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two control volumes
    ///
    /// the thermal conductance is determined by 
    /// Thermal conductance 
    /// /// (2 * pi * L * K)/
    /// ln(outer_radius/inner_radius)
    ///
    /// using obtain_thermal_conductance_annular_cylinder
    /// under common_functions
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
    CylindricalConductionConvectionLiquidOutside(
        (Material,RadialCylindricalThicknessThermalConduction,
        ThermodynamicTemperature,Pressure),
        (HeatTransfer, 
        OuterDiameterThermalConduction, 
        CylinderLengthThermalConduction),
    ),

    /// 1D Cylindrical Coordinates Thermal Resistance
    /// We return a ThermalConductance because it's more convenient
    ///
    /// basically have three control volumes along the outer wall
    ///
    /// -------------------------------------------------------> r
    /// //                           ----------------------------
    /// //                           |                          |                          
    /// // *                         *         solid_cv_1       *                  
    /// //                           |                          |                   
    /// // fluid node                ----------------------------
    /// // (T_f)                solid_surface
    ///
    /// Where r is the radius 
    /// basically the liquid is on the inside (larger smaller r)
    ///
    /// between solid_cv_1 and solid_surface 
    /// there is a thermal resistance 
    /// based on a q'' = k dT/dr
    ///
    /// between solid_surface and fluid_node, there is convection resistance
    /// specified by a Nusselt Number so that we get a heat transfer 
    /// coefficient
    ///
    /// For the conduction bit,
    /// we have one material which determines conductivity 
    /// and then length which determines the distance between 
    /// the two control volumes
    ///
    /// the thermal conductance is determined by 
    /// Thermal conductance 
    /// /// (2 * pi * L * K)/
    /// ln(outer_radius/inner_radius)
    ///
    /// using obtain_thermal_conductance_annular_cylinder
    /// under common_functions
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
    CylindricalConductionConvectionLiquidInside(
        (Material,RadialCylindricalThicknessThermalConduction,
        ThermodynamicTemperature,Pressure),
        (HeatTransfer, 
        InnerDiameterThermalConduction, 
        CylinderLengthThermalConduction),
    ),


    /// The user Specifies a heat Addition for the BC
    /// The uom type is Power
    UserSpecifiedHeatAddition,
    
    /// Use this enum to specify a constant heat flux
    /// you will, of course, need to provide an area
    UserSpecifiedHeatFluxCustomArea(Area),

    /// Use this enum to identify that you are 
    /// specifying a curved cylindrical surface area 
    /// on the outer surface of a cylinder
    UserSpecifiedHeatFluxCylindricalOuterArea(
        CylinderLengthThermalConduction,
        OuterDiameterThermalConduction,
    ),

    /// Use this enum to identify that you are 
    /// specifying a curved cylindrical surface area 
    /// on the inner surface of a cylinder
    UserSpecifiedHeatFluxCylindricalInnerArea(
        CylinderLengthThermalConduction,
        InnerDiameterThermalConduction,
    ),

    
    /// For convection between solid and fluid, 
    /// the heat flux from solid surface to fluid 
    /// is:
    ///
    /// q = h A(T_s - T_f)
    /// 
    /// this interaction calculates power based on a given h and A 
    UserSpecifiedConvectionResistance(
        DataUserSpecifiedConvectionResistance),


    /// For advection one would only specify the mass flowrate 
    /// from one control volume to another
    Advection(DataAdvection),
}

/// here we have a struct for simple convection resistance
/// in three dimensions
/// on
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataUserSpecifiedConvectionResistance{

    /// surface area for heat convection
    pub surf_area: SurfaceArea,
    /// heat transfer coefficient in watts per meter per kelvin
    pub heat_transfer_coeff: HeatTransfer,

}

/// here we have a useful for necessary advection information 

#[derive(Debug,Clone,Copy,PartialEq)]
pub struct DataAdvection{

    /// mass flowrate
    pub mass_flowrate: MassRate,
    /// fluid density of control volume on left
    ///
    /// which means when you link control volumes or boundary 
    /// link(cv1, cv2, interaction)
    ///
    /// the picture is like this 
    ///
    /// (cv1) ----> advection ---> (cv2)
    ///
    /// cv1 is the left control volume 
    /// cv2 is the right control volume
    ///
    /// now, the cv is not always a cv, it could be any heat
    /// transfer entity
    pub fluid_density_heat_transfer_entity_1: MassDensity,
    /// fluid density of control volume on left
    ///
    /// which means when you link control volumes or boundary 
    /// link(cv1, cv2, interaction)
    ///
    /// the picture is like this 
    ///
    /// (cv1) ----> advection ---> (cv2)
    ///
    /// cv1 is the left control volume 
    /// cv2 is the right control volume
    /// now, the cv is not always a cv, it could be any heat
    /// transfer entity
    pub fluid_density_heat_transfer_entity_2: MassDensity

}


impl Into<HeatTransferInteractionType> for DataAdvection {
    fn into(self) -> HeatTransferInteractionType {
        HeatTransferInteractionType::Advection(self)
    }
}
