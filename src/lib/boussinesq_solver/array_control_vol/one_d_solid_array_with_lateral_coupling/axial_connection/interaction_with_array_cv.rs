use crate::boussinesq_solver::array_control_vol::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::boussinesq_solver::single_control_vol::SingleCVNode;
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
}
