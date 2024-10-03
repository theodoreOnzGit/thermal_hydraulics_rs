use crate::tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::tuas_boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::tuas_boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::tuas_boussinesq_solver::single_control_vol::SingleCVNode;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

impl SolidColumn {

    /// attaches an array control volume to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    ///
    /// in this case, it is a solid column
    pub fn link_solid_column_to_the_front_of_this_solid_column(
        &mut self,
        solid_column_other: &mut SolidColumn,
        interaction: HeatTransferInteractionType,) -> Result<(), ThermalHydraulicsLibError>{

        if let HeatTransferInteractionType::Advection(_) = interaction {
            println!("You cannot have advection interactions for Solid Columns");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
        }
        // basically we need to get the front of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.front_single_cv;

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        &mut solid_column_other.back_single_cv;

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_self,
            single_cv_node_other,
            interaction)
    }

    /// attaches an array control volume to the back of this 
    /// array control volume 
    /// (back --- cv_other --- front) ---- (back --- cv_self --- front)
    ///
    /// in this case, it is a solid column
    pub fn link_solid_column_to_the_back_of_this_solid_column(
        &mut self,
        solid_column_other: &mut SolidColumn,
        interaction: HeatTransferInteractionType,) -> Result<(), ThermalHydraulicsLibError>{

        if let HeatTransferInteractionType::Advection(_) = interaction {
            println!("You cannot have advection interactions for Solid Columns");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
        }
        // basically we need to get the back of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.back_single_cv;

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        &mut solid_column_other.front_single_cv;

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches an solid column array control volume to the front of this 
    /// fluid array control volume 
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    ///
    pub fn link_fluid_array_to_the_front_of_this_solid_column(
        &mut self,
        fluid_array_other: &mut FluidArray,
        interaction: HeatTransferInteractionType,) -> Result<(), ThermalHydraulicsLibError>{

        if let HeatTransferInteractionType::Advection(_) = interaction {
            println!("You cannot have advection interactions for Solid Columns");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
        }
        // basically we need to get the front of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.front_single_cv;

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        &mut fluid_array_other.back_single_cv;

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_self,
            single_cv_node_other,
            interaction)
    }

    /// attaches an solid column 
    /// array control volume to the back of this 
    /// fluid array control volume 
    /// (back --- cv_other --- front) ---- (back --- cv_self --- front)
    ///
    pub fn link_fluid_array_to_the_back_of_this_solid_column(
        &mut self,
        fluid_array_other: &mut FluidArray,
        interaction: HeatTransferInteractionType,) -> Result<(), ThermalHydraulicsLibError>{

        if let HeatTransferInteractionType::Advection(_) = interaction {
            println!("You cannot have advection interactions for Solid Columns");
            return Err(ThermalHydraulicsLibError::WrongHeatTransferInteractionType);
        }
        // basically we need to get the back of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.back_single_cv;

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        &mut fluid_array_other.front_single_cv;

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }
}
