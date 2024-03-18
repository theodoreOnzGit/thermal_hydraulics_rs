use crate::boussinesq_solver::array_control_vol::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::single_control_vol::SingleCVNode;
use crate::boussinesq_solver::heat_transfer_correlations::heat_transfer_interactions::heat_transfer_interaction_enums::HeatTransferInteractionType;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;

use uom::si::f64::*;

impl FluidArray {

    /// attaches a single cv to the front,entrance,
    /// lower or inner side of the 
    /// array cv
    /// 
    ///
    /// basically in whatever coordinate system, it is the lowest 
    /// value 
    ///
    /// for spheres, the lowest r (inner side)
    /// for cylinders, the lowest r or z (inner or lower)
    /// for cartesian, the lowest x, y or z (back)
    ///
    /// We use this convention because heat can flow either way 
    /// and fluid can flow either way. 
    ///
    /// So it's better to define higher and lower based upon a coordinate 
    /// axis
    pub fn link_single_cv_to_lower_side(&mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError>{


        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = &mut self.back_single_cv;

        // now link both cvs or calculate between them

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches a single cv to the exit,back,
    /// higher or outer side of the 
    /// array cv
    /// 
    ///
    /// basically in whatever coordinate system, it is the lowest 
    /// value 
    ///
    /// for spheres, the highest r (outer side)
    /// for cylinders, the highest r or z (outer or higher)
    /// for cartesian, the highest x, y or z (front)
    ///
    /// We use this convention because heat can flow either way 
    /// and fluid can flow either way. 
    ///
    /// So it's better to define higher and lower based upon a coordinate 
    /// axis
    pub fn link_single_cv_to_higher_side(&mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<(), ThermalHydraulicsLibError>{

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = &mut self.front_single_cv;

        // now link both cvs or calculate between them

        SingleCVNode::calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// calculates timestep for a single cv attached to the front of the 
    /// array cv
    /// (back --- cv_self --- front) ---- (single cv)
    pub fn calculate_timestep_for_single_cv_to_front_of_array_cv(
        &mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.front_single_cv;

        single_cv_node_self.calculate_mesh_stability_timestep_for_two_single_cv_nodes(
            single_cv_node_other,
            interaction)

    }

    /// calculates timestep for a single cv attached to the back of the 
    /// array cv
    /// (single cv) --- (back --- cv_self --- front)
    pub fn calculate_timestep_for_single_cv_to_back_of_array_cv(
        &mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,ThermalHydraulicsLibError> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        &mut self.back_single_cv;

        single_cv_node_self.calculate_mesh_stability_timestep_for_two_single_cv_nodes(
            single_cv_node_other,
            interaction)
    }
}
