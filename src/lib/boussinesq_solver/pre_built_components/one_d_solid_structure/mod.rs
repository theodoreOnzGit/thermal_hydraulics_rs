use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;

use super::heat_transfer_entities::cv_types::CVType;
use super::heat_transfer_entities::HeatTransferEntity;
use uom::si::f64::*;

/// The simplest component is an insulated pipe
///
/// This is a simple pipe with a set hydraulic diameter and length
///
/// the standard assumption is that at each boundary of this pipe,
/// there is no conduction heat transfer in the axial direction
#[derive(Clone,Debug,PartialEq)]
pub struct SolidStructure {

    inner_nodes: usize,

    /// this HeatTransferEntity represents the pipe shell 
    /// only one radial layer of control volumes is used to simulate 
    /// the pipe shell
    ///
    /// it is thermally coupled to insulation and to the fluid 
    /// in the pipe_fluid_array
    pub solid_array: HeatTransferEntity,


    /// this HeatTransferEntity represents the pipe fluid
    /// which is coupled to the pipe shell via a Nusselt Number based
    /// thermal resistance (usually Gnielinski correlation)
    pub pipe_fluid_array: HeatTransferEntity,

    /// this HeatTransferEntity represents the pipe insulation
    ///
    /// which is 
    /// exposed to an ambient constant temperature boundary condition
    /// This is because constant heat flux BCs are not common for pipes
    /// except for fully/ideally insulated pipes
    ///
    /// this is coupled to the pipe_shell
    pub insulation: HeatTransferEntity,

    /// pipe ambient temperature
    pub ambient_temperature: ThermodynamicTemperature,

    /// pipe heat transfer coefficient to ambient
    pub heat_transfer_to_ambient: HeatTransfer,

    /// pipe outer diameter (tube)
    pub tube_od: Length,

    /// pipe inner diameter (tube)
    pub tube_id: Length,

    /// pipe outer diameter (insulation)
    insulation_od: Length,

    /// pipe inner diameter (insulation)
    insulation_id: Length,

    /// flow area
    flow_area: Area,

    /// loss correlations
    darcy_loss_correlation: DimensionlessDarcyLossCorrelations,

}

impl SolidStructure {

    /// constructs a new insulated pipe
    ///
    /// you need to supply the initial temperature, ambient temperature
    /// as well as all the pipe parameters 
    ///
    /// such as:
    ///
    /// 1. flow area 
    /// 2. hydraulic diameter 
    /// 3. incline angle
    /// 4. any form losses beyond the Gnielinski correlation
    /// 5. inner diameter (id)
    /// 6. shell outer diameter (od) assumed to be same as insulation id
    /// 7. pipe shell material 
    /// 8. pipe fluid 
    /// 9. fluid pressure (if in doubt, 1 atmosphere will do)
    /// 10. solid pressure (if in doubt, 1 atmosphere will do)
    /// 11. heat transfer coeffficient to ambient
    /// 12. how many inner axial nodes for both solid and fluid arrays
    /// 13. insulation thickness 
    /// 14. darcy loss correlation
    ///
    /// The number of total axial nodes is the number of inner nodes plus 2
    ///
    /// this is because there are two nodes at the periphery of the pipe 
    /// and there
    /// at each timestep, you are allowed to set a heater power, where 
    /// heat is dumped into the heated tube surrounding the pipe
    ///
    /// so the pipe shell becomes the heating element so to speak
    pub fn new_insulated_pipe(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature,
        fluid_pressure: Pressure,
        solid_pressure: Pressure,
        flow_area: Area,
        incline_angle: Angle,
        form_loss: Ratio,
        shell_id: Length,
        shell_od: Length,
        insulation_thickness: Length,
        pipe_length: Length,
        hydraulic_diameter: Length,
        pipe_shell_material: SolidMaterial,
        pipe_fluid: LiquidMaterial,
        htc_to_ambient: HeatTransfer,
        user_specified_inner_nodes: usize,
        surface_roughness: Length) -> SolidStructure {

        // inner fluid_array
        let fluid_array: FluidArray = 
        FluidArray::new_odd_shaped_pipe(
            pipe_length,
            hydraulic_diameter,
            flow_area,
            initial_temperature,
            fluid_pressure,
            pipe_shell_material,
            pipe_fluid,
            form_loss,
            user_specified_inner_nodes,
            incline_angle
        );

        // now the outer pipe array
        let pipe_shell = 
        SolidColumn::new_cylindrical_shell(
            pipe_length,
            shell_id,
            shell_od,
            initial_temperature,
            solid_pressure,
            pipe_shell_material,
            user_specified_inner_nodes 
        );

        let insulation_id = shell_od;
        let insulation_od = insulation_id + insulation_thickness;

        // insulation
        let insulation = 
        SolidColumn::new_cylindrical_shell(
            pipe_length,
            insulation_id,
            insulation_od,
            initial_temperature,
            solid_pressure,
            SolidMaterial::Fiberglass,
            user_specified_inner_nodes 
        );

        // fluid pipe loss correlation
        //


        let pipe_loss_correlation = DimensionlessDarcyLossCorrelations::
            new_pipe(pipe_length, surface_roughness, hydraulic_diameter, form_loss);

        return Self { inner_nodes: user_specified_inner_nodes,
            solid_array: CVType::SolidArrayCV(pipe_shell).into(),
            pipe_fluid_array: CVType::FluidArrayCV(fluid_array).into(),
            ambient_temperature,
            heat_transfer_to_ambient: htc_to_ambient,
            tube_od: shell_od,
            tube_id: shell_id,
            insulation_od: shell_od,
            insulation_id: shell_od+insulation_thickness,
            flow_area,
            darcy_loss_correlation: pipe_loss_correlation,
            insulation: CVType::SolidArrayCV(insulation).into(),
        };
    }
}


/// stuff such as conductances are calculated here
pub mod preprocessing;

/// implementations for the FluidComponent trait
/// are done here
pub mod fluid_component;


/// stuff for calculation is done here, ie, advancing timestep
pub mod calculation;

/// postprocessing stuff, ie, get the temperature vectors 
/// of both arrays of control volumes 
pub mod postprocessing;
