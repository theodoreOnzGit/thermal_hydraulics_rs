
use uom::si::f64::*;
use crate::heat_transfer_lib::thermophysical_properties::Material;
use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::*;


use crate::heat_transfer_lib::control_volume_calculations:: 
heat_transfer_entities::HeatTransferEntity;
/// basically an enum for you to specify 
/// if the liquid on the inner curved surface of the shell or outer 
/// curved surface of the shell
///
/// in the context of a convection and conductivity 
/// thermal resistance calculation,
///
#[derive(Debug,Clone,Copy,PartialEq)]
pub enum CylindricalAndSphericalSolidFluidArrangement {
    /// indicates that fluid in the inner side of a curved shell 
    ///
    /// -----------------------------------------> r
    /// fluid               ||                  solid
    ///
    FluidOnInnerSurfaceOfSolidShell,
    /// indicates that fluid in the inner side of a curved shell 
    ///
    /// -----------------------------------------> r
    /// solid               ||                  fluid
    ///
    FluidOnOuterSurfaceOfSolidShell
}

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
    UserSpecifiedHeatAddition(Power),
}
