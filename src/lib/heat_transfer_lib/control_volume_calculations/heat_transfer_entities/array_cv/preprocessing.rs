use std::f64::consts::PI;

use crate::fluid_mechanics_lib::prelude::FluidComponent;
use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::ArrayCVType;
use crate::heat_transfer_lib::
control_volume_calculations::heat_transfer_entities::SingleCVNode;
use uom::si::f64::*;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::enum_selection_alpha::
interactions_single_cv::*;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::HeatTransferInteractionType;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::enum_selection_alpha::
interactions_single_cv::calculate_between_two_singular_cv_nodes;

use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::
enum_selection_alpha::
timestep_control_single_cv::calculate_mesh_stability_conduction_timestep_for_single_node_and_bc;
use crate::heat_transfer_lib::control_volume_calculations::
heat_transfer_interactions::
enum_selection_alpha::
timestep_control_single_cv::calculate_mesh_stability_timestep_for_two_single_cv_nodes;


use crate::heat_transfer_lib::thermophysical_properties::
Material::{Solid, Liquid};



impl ArrayCVType {
    /// gets the maximum timestep for the arrayCV
    pub fn get_max_timestep(
        &mut self,
        max_temperature_change: TemperatureInterval) 
    -> Result<Time, String>{
        match self {
            ArrayCVType::Cartesian1D(cv) => {
                cv.get_max_timestep(max_temperature_change)
            },
            ArrayCVType::GenericPipe(cv) => {
                let mass_flowrate = cv.get_mass_flowrate();

                Ok(cv.get_max_timestep(
                max_temperature_change,
                mass_flowrate).unwrap())
            },

        }
    }
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
        interaction: HeatTransferInteractionType) -> Result<(), String>{


        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_back_cv_deref_mut().unwrap();

        // now link both cvs or calculate between them

        calculate_between_two_singular_cv_nodes(
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
        interaction: HeatTransferInteractionType) -> Result<(), String>{

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();

        // now link both cvs or calculate between them

        calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches an array control volume to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    pub fn link_array_cv_to_the_front_of_this_cv(
        &mut self,
        array_cv_other: &mut ArrayCVType,
        interaction: HeatTransferInteractionType,) -> Result<(), String>{

        // basically we need to get the front of the self cv, 
        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();

        // and the back of the other cv
        let single_cv_node_other: &mut SingleCVNode = 
        array_cv_other.nested_back_cv_deref_mut().unwrap();

        calculate_between_two_singular_cv_nodes(
            single_cv_node_other,
            single_cv_node_self,
            interaction)
    }

    /// attaches an constant heat flux BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant heat flux)
    pub fn link_heat_flux_bc_to_front_of_this_cv(
        &mut self,
        heat_flux_into_control_vol: HeatFluxDensity,
        interaction: HeatTransferInteractionType) -> Result<(),String>{

        // first, obtain a heat transfer area from the constant heat flux 
        // BC
        let heat_transfer_area: Area = match interaction{
            HeatTransferInteractionType::
                UserSpecifiedThermalConductance(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                SingleCartesianThermalConductanceOneDimension(_, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                DualCartesianThermalConductance(_, _) => 
                return Err("please specify interaction \n 
                    type as UserSpecifiedHeatFluxCustomArea \n 
                    or Similar".to_string()),

            HeatTransferInteractionType::
                DualCylindricalThermalConductance(_, _, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or \n
                    Similar".to_string()),

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidOutside(_, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(_, _) => 
                return Err("please specify interaction \n 
                    type as UserSpecifiedHeatFluxCustomArea \n 
                    or Similar".to_string()),

            HeatTransferInteractionType::
                UserSpecifiedHeatAddition => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                DualCartesianThermalConductanceThreeDimension(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
            // these interaction types are acceptable
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCustomArea(area) => area,

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(
                cylinder_length, od) => {
                    let od: Length = od.into();
                    let cylinder_length: Length  = cylinder_length.into();

                    let area = PI * od * cylinder_length;
                    area
                },

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalInnerArea(
                cylinder_length, id) => {
                    let id: Length = id.into();
                    let cylinder_length: Length  = cylinder_length.into();

                    let area = PI * id * cylinder_length;
                    area

                },

            HeatTransferInteractionType::
                UserSpecifiedConvectionResistance(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
            HeatTransferInteractionType::
                Advection(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        };

        // once area is calculated, we can calculate heat flowrate into 
        // cv
        let heat_flowrate_into_control_vol: Power = 
        heat_flux_into_control_vol * heat_transfer_area;

        // then push the power to the front control volume, 
        // then just do conduction timescales 
        // There's no advection in this case so no need to worry 

        let front_cv_ref = self.nested_front_cv_deref_mut().unwrap();

        // push the power power to the front cv 

        front_cv_ref.rate_enthalpy_change_vector.push( 
            heat_flowrate_into_control_vol);
        // calculate conduction timescales
        calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
            front_cv_ref,
            interaction)?;

        // I don't calculate solid-liquid timescales here, 
        // kind of redundant 
        // I could implement it in future though
        // todo: nusselt number adjusted timescales for heatflux 


        return Ok(());
    }

    /// attaches an constant heat flux BC to the front of this 
    /// array control volume 
    /// (constant heat flux) ---- (back --- cv_self --- front)
    pub fn link_heat_flux_bc_to_back_of_this_cv(
        &mut self,
        heat_flux_into_control_vol: HeatFluxDensity,
        interaction: HeatTransferInteractionType) -> Result<(),String>{

        // first, obtain a heat transfer area from the constant heat flux 
        // BC
        let heat_transfer_area: Area = match interaction{
            HeatTransferInteractionType::
                UserSpecifiedThermalConductance(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                SingleCartesianThermalConductanceOneDimension(_, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                DualCartesianThermalConductance(_, _) => 
                return Err("please specify interaction \n 
                    type as UserSpecifiedHeatFluxCustomArea \n 
                    or Similar".to_string()),

            HeatTransferInteractionType::
                DualCylindricalThermalConductance(_, _, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or \n
                    Similar".to_string()),

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidOutside(_, _) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                CylindricalConductionConvectionLiquidInside(_, _) => 
                return Err("please specify interaction \n 
                    type as UserSpecifiedHeatFluxCustomArea \n 
                    or Similar".to_string()),

            HeatTransferInteractionType::
                UserSpecifiedHeatAddition => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),

            HeatTransferInteractionType::
                DualCartesianThermalConductanceThreeDimension(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
            // these interaction types are acceptable
            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCustomArea(area) => area,

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalOuterArea(
                cylinder_length, od) => {
                    let od: Length = od.into();
                    let cylinder_length: Length  = cylinder_length.into();

                    let area = PI * od * cylinder_length;
                    area
                },

            HeatTransferInteractionType::
                UserSpecifiedHeatFluxCylindricalInnerArea(
                cylinder_length, id) => {
                    let id: Length = id.into();
                    let cylinder_length: Length  = cylinder_length.into();

                    let area = PI * id * cylinder_length;
                    area

                },

            HeatTransferInteractionType::
                UserSpecifiedConvectionResistance(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
            HeatTransferInteractionType::
                Advection(_) => 
                return Err("please specify interaction type as \n 
                    UserSpecifiedHeatFluxCustomArea or Similar".to_string()),
        };

        // once area is calculated, we can calculate heat flowrate into 
        // cv
        let heat_flowrate_into_control_vol: Power = 
        heat_flux_into_control_vol * heat_transfer_area;

        // then push the power to the front control volume, 

        // then push the power to the front control volume, 
        // then just do conduction timescales 
        // There's no advection in this case so no need to worry 

        let back_cv_ref = self.nested_front_cv_deref_mut().unwrap();

        // push the power power to the back cv 

        back_cv_ref.rate_enthalpy_change_vector.push( 
            heat_flowrate_into_control_vol);
        // calculate conduction timescales
        calculate_mesh_stability_conduction_timestep_for_single_node_and_bc(
            back_cv_ref,
            interaction)?;

        // I don't calculate solid-liquid timescales here, 
        // kind of redundant 
        // I could implement it in future though
        // todo: nusselt number adjusted timescales for heatflux 

        return Ok(());
    }

    /// attaches an constant heat rate BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant heat rate)
    pub fn link_heat_addition_to_front_of_this_cv(
        &mut self,
        heat_rate: Power,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();

        // then link it to the other bc as per normal for 
        // single control volumes

        calculate_constant_heat_addition_front_single_cv_back(
            single_cv_node_self,
            heat_rate,
            interaction
        )

    }

    /// attaches an constant heat rate BC to the front of this 
    /// array control volume 
    /// (constant heat rate) --- (back --- cv_self --- front)
    pub fn link_heat_addition_to_back_of_this_cv(
        &mut self,
        heat_rate: Power,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_back_cv_deref_mut().unwrap();

        calculate_single_cv_front_constant_heat_addition_back(
            heat_rate,
            single_cv_node_self,
            interaction
        )
    }

    /// attaches an constant temperature BC to the front of this 
    /// array control volume 
    /// (back --- cv_self --- front) ---- (constant temperature)
    pub fn link_constant_temperature_to_front_of_this_cv(
        &mut self,
        bc_temperature: ThermodynamicTemperature,
        interaction: HeatTransferInteractionType) -> Result<(),String> {
        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();


        calculate_constant_temperature_front_single_cv_back (
            single_cv_node_self,
            bc_temperature,
            interaction
        )
    }

    /// attaches an constant temperature BC to the front of this 
    /// array control volume 
    /// (constant temperature) --- (back --- cv_self --- front)
    pub fn link_constant_temperature_to_back_of_this_cv(
        &mut self,
        bc_temperature: ThermodynamicTemperature,
        interaction: HeatTransferInteractionType) -> Result<(),String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_back_cv_deref_mut().unwrap();

        calculate_single_cv_node_front_constant_temperature_back (
            bc_temperature,
            single_cv_node_self,
            interaction
        )
    }
    /// calculates timestep for a single cv attached to the front of the 
    /// array cv
    /// (back --- cv_self --- front) ---- (single cv)
    pub fn calculate_timestep_for_single_cv_to_front_of_array_cv(
        &mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();

        calculate_mesh_stability_timestep_for_two_single_cv_nodes(
            single_cv_node_self,
            single_cv_node_other,
            interaction)

    }

    /// calculates timestep for a single cv attached to the back of the 
    /// array cv
    /// (single cv) --- (back --- cv_self --- front)
    pub fn calculate_timestep_for_single_cv_to_back_of_array_cv(
        &mut self,
        single_cv_node_other: &mut SingleCVNode,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_back_cv_deref_mut().unwrap();

        calculate_mesh_stability_timestep_for_two_single_cv_nodes(
            single_cv_node_self,
            single_cv_node_other,
            interaction)
    }


    /// calculates timestep for an array cv attached to the front of the 
    /// array cv
    /// (back --- cv_self --- front) ---- (back --- cv_other --- front)
    pub fn calculate_timestep_for_array_cv_to_front_of_this_array_cv(
        &mut self,
        array_cv_other: &mut ArrayCVType,
        interaction: HeatTransferInteractionType) -> Result<Time,String> {

        // we need to obtain the single cv from the array cv first 
        // and this will be the front cv or outer cv 
        //

        let single_cv_node_self: &mut SingleCVNode = 
        self.nested_front_cv_deref_mut().unwrap();

        // we need to obtain the single cv from the array cv first 
        // and this will be the back cv or inner cv 
        //

        let single_cv_node_other: &mut SingleCVNode = 
        array_cv_other.nested_back_cv_deref_mut().unwrap();

        calculate_mesh_stability_timestep_for_two_single_cv_nodes(
            single_cv_node_self,
            single_cv_node_other,
            interaction)
    }
}
