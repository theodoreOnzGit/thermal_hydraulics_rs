use std::f64::consts::PI;

use super::SimpleShellAndTubeHeatExchanger;
use uom::si::f64::*;
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
use uom::si::thermodynamic_temperature::kelvin;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::boussinesq_solver::heat_transfer_correlations::nusselt_number_correlations::enums::NusseltCorrelation;
use crate::boussinesq_solver::heat_transfer_correlations::nusselt_number_correlations::input_structs::GnielinskiData;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::SolidMaterial;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::LiquidMaterial;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;

impl SimpleShellAndTubeHeatExchanger {

    /// gets the temperature of the inner pipe shell array
    pub fn inner_pipe_shell_temperature(&mut self) -> 
        Result<Vec<ThermodynamicTemperature>, ThermalHydraulicsLibError>{
            self.inner_pipe_shell_array_for_single_tube.get_temperature_vector()
    }

    /// gets the temperature of the tube side fluid array
    pub fn inner_tube_fluid_array_temperature(&mut self) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.inner_pipe_shell_array_for_single_tube.get_temperature_vector()
    }

    /// gets the shell side fluid array temperature
    pub fn shell_side_fluid_array_temperature(&mut self,) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.shell_side_fluid_array.get_temperature_vector()
    }

    /// gets the shell side outer tube temperature 
    pub fn shell_side_outer_tube_array_temperature(&mut self,) -> 
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
            self.outer_shell.get_temperature_vector()
    }

    /// gets the temperature of the insulation 
    pub fn insulation_array_temperature(&mut self,) ->
        Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{

            if self.heat_exchanger_has_insulation {
                self.insulation_array.get_temperature_vector()
            } else {
                // if the heat exchanger is not even insulated, then 
                // there should be no insulation array
                todo!("simple Shell and Tube Heat exchanger is not insulated")
            }
    }

    /// provides the overall heat transfer coefficient based on the 
    /// shell side area
    ///
    /// note that this assumes the tubes are circular, not oddly 
    /// shaped
    ///
    /// as mentioned by Du:
    ///
    /// 1/U = 1/h_t d_o/d_i + d_o/(2 lambda_w) ln (d_o/d_i) 
    /// + 1/h_s
    pub fn overall_heat_transfer_coeff_u_shell_side(&self,
        correct_for_prandtl_wall_temperatures: bool) 
        -> Result<HeatTransfer, ThermalHydraulicsLibError>{

            // the thermal conductance here should be based on the 
            // nusselt number correlation

            // before any calculations, I will first need a clone of 
            // the fluid array and inner shell array
            //
            // the fluid array represents only a single tube
            let mut tube_side_single_fluid_array_clone: FluidArray = 
                self.tube_side_fluid_array_for_single_tube.clone().try_into()?;


            let mut pipe_shell_clone: SolidColumn = 
                self.inner_pipe_shell_array_for_single_tube.clone().try_into()?;

            // also need to get basic temperatures and mass flowrates 
            // only do this once because some of these methods involve 
            // cloning, which is computationally expensive

            let single_tube_mass_flowrate: MassRate = 
                tube_side_single_fluid_array_clone.get_mass_flowrate();

            let tube_fluid_temperature: ThermodynamicTemperature 
                = tube_side_single_fluid_array_clone.try_get_bulk_temperature()?;

            let wall_temperature: ThermodynamicTemperature 
                = pipe_shell_clone.try_get_bulk_temperature()?;

            let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

            let pipe_shell_surf_temperature: ThermodynamicTemperature 
                = pipe_shell_clone.try_get_bulk_temperature()?;

            let single_tube_hydraulic_diameter = 
                self.get_tube_side_hydraulic_diameter();
            let single_tube_flow_area: Area = 
                tube_side_single_fluid_array_clone.get_cross_sectional_area_immutable();

            // flow area and hydraulic diameter are ok


            let tube_fluid_material: LiquidMaterial
                = tube_side_single_fluid_array_clone.material_control_volume.try_into()?;

            let solid_material: SolidMaterial 
                = pipe_shell_clone.material_control_volume.try_into()?;

            let viscosity: DynamicViscosity = 
                tube_fluid_material.try_get_dynamic_viscosity(tube_fluid_temperature)?;

            // need to convert hydraulic diameter to an equivalent 
            // spherical diameter
            //
            // but for now, I'm going to use Re and Nu using hydraulic diameter 
            // and live with it for the time being
            //
            let reynolds_number_single_tube_abs_for_nusselt_estimate: Ratio = 
                (single_tube_mass_flowrate/
                single_tube_flow_area
                *single_tube_hydraulic_diameter / viscosity).abs();

            // next, bulk prandtl number 

            let bulk_prandtl_number: Ratio 
                = tube_fluid_material.try_get_prandtl_liquid(
                    tube_fluid_temperature,
                    atmospheric_pressure
                )?;

            let nusselt_estimate_tube_side: Ratio; 
            let darcy_friction_factor_tube_side: Ratio = self.
                tube_side_custom_component_loss_correlation.
                darcy_friction_factor_fldk(reynolds_number_single_tube_abs_for_nusselt_estimate)
                .unwrap();

            if correct_for_prandtl_wall_temperatures {

                // then wall prandtl number (partially corrected)

                let part_correct_wall_temperature: ThermodynamicTemperature = 
                    ThermodynamicTemperature::new::<kelvin>(
                        0.1 * (
                            3.0 * wall_temperature.get::<kelvin>() + 
                            7.0 * tube_fluid_temperature.get::<kelvin>()
                        )
                    );

                let wall_prandtl_number_part_correct: Ratio 
                    = tube_fluid_material.try_get_prandtl_liquid(
                        part_correct_wall_temperature,
                        atmospheric_pressure
                    )?;

                nusselt_estimate_tube_side = 
                    self.tube_side_nusselt_correlation
                    .estimate_based_on_prandtl_darcy_and_reynolds_wall_correction(
                        bulk_prandtl_number, 
                        wall_prandtl_number_part_correct, 
                        darcy_friction_factor_tube_side,
                        reynolds_number_single_tube_abs_for_nusselt_estimate)?;
            } else {
                nusselt_estimate_tube_side = 
                    self.tube_side_nusselt_correlation
                    .estimate_based_on_prandtl_darcy_and_reynolds_wall_correction(
                        bulk_prandtl_number, 
                        bulk_prandtl_number, 
                        darcy_friction_factor_tube_side,
                        reynolds_number_single_tube_abs_for_nusselt_estimate)?;
            }

            let k_fluid_average: ThermalConductivity = 
                tube_fluid_material.try_get_thermal_conductivity(
                    tube_fluid_temperature)?;

            let tube_side_h
                = nusselt_estimate_tube_side 
                * k_fluid_average / single_tube_hydraulic_diameter;
            
            let d_o_by_d_i: Ratio = 
                self.tube_side_od/self.tube_side_id;

            // 1/h_t * d_o/d_i 
            let one_over_ht_times_do_by_di = 
                1.0 as f64 / tube_side_h 
                * d_o_by_d_i;
            
            // now for the shell conductivity part 
            let shell_conductivity_average: ThermalConductivity = 
                solid_material.try_get_thermal_conductivity(
                    pipe_shell_surf_temperature)?;

            let ln_d_o_by_d_i: f64 = 
                d_o_by_d_i.get::<ratio>().ln();

            // d_o/(2 lambda_w) * ln (d_o_by_d_i) 

            let do_by_2_lambda_w_times_ln_do_by_di 
                = self.tube_side_od * 0.5 /
                (shell_conductivity_average) * 
                ln_d_o_by_d_i;

            // now for shell side heat transfer coeff

            // the thermal conductance here should be based on the 
            // nusselt number correlation

            // before any calculations, I will first need a clone of 
            // the fluid array and twisted tape array
            let mut shell_side_fluid_array_clone: FluidArray = 
                self.shell_side_fluid_array.clone().try_into()?;

            // also need to get basic temperatures and mass flowrates 
            // only do this once because some of these methods involve 
            // cloning, which is computationally expensive

            let shell_side_mass_flowrate: MassRate = 
                shell_side_fluid_array_clone.get_mass_flowrate();

            let shell_side_fluid_temperature: ThermodynamicTemperature 
                = shell_side_fluid_array_clone.try_get_bulk_temperature()?;

            let wall_temperature: ThermodynamicTemperature 
                = pipe_shell_clone.try_get_bulk_temperature()?;

            let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);

            let shell_side_fluid_hydraulic_diameter = 
                self.get_shell_side_hydraulic_diameter();

            let shell_side_cross_sectional_flow_area: Area = 
                self.get_shell_side_cross_sectional_area();


            // flow area and hydraulic diameter are ok


            let fluid_material: LiquidMaterial
                = shell_side_fluid_array_clone.material_control_volume.try_into()?;

            let viscosity: DynamicViscosity = 
                fluid_material.try_get_dynamic_viscosity(shell_side_fluid_temperature)?;

            // need to convert hydraulic diameter to an equivalent 
            // spherical diameter
            //
            // but for now, I'm going to use Re and Nu using hydraulic diameter 
            // and live with it for the time being
            //
            let reynolds_number_shell_side_abs_for_nusselt_estimate: Ratio = 
                (shell_side_mass_flowrate/
                shell_side_cross_sectional_flow_area
                *shell_side_fluid_hydraulic_diameter / viscosity).abs();

            // next, bulk prandtl number 

            let bulk_prandtl_number: Ratio 
                = fluid_material.try_get_prandtl_liquid(
                    shell_side_fluid_temperature,
                    atmospheric_pressure
                )?;



            let shell_side_fluid_to_inner_tube_surf_nusselt_correlation: NusseltCorrelation
                = self.shell_side_nusselt_correlation_to_tubes;
            let darcy_friction_factor_shell_side: Ratio = self.
                shell_side_custom_component_loss_correlation.
                darcy_friction_factor_fldk(reynolds_number_shell_side_abs_for_nusselt_estimate)
                .unwrap();

            let mut pipe_prandtl_reynolds_gnielinksi_data: GnielinskiData 
                = GnielinskiData::default();
            pipe_prandtl_reynolds_gnielinksi_data.reynolds = reynolds_number_shell_side_abs_for_nusselt_estimate;
            pipe_prandtl_reynolds_gnielinksi_data.prandtl_bulk = bulk_prandtl_number;
            pipe_prandtl_reynolds_gnielinksi_data.prandtl_wall = bulk_prandtl_number;
            pipe_prandtl_reynolds_gnielinksi_data.darcy_friction_factor = 
                darcy_friction_factor_shell_side;
            pipe_prandtl_reynolds_gnielinksi_data.length_to_diameter = 
                shell_side_fluid_array_clone.get_component_length_immutable()/
                shell_side_fluid_hydraulic_diameter;


            // I need to use Nusselt correlations present in this struct 
            //
            // wall correction is optionally done here
            //
            // this uses the gnielinski correlation for pipes or tubes

            let nusselt_estimate_shell_side: Ratio;

            if correct_for_prandtl_wall_temperatures {

                // then wall prandtl number (partially corrected)

                let part_correct_wall_temperature: ThermodynamicTemperature = 
                    ThermodynamicTemperature::new::<kelvin>(
                        0.1 * (
                            3.0 * wall_temperature.get::<kelvin>() + 
                            7.0 * shell_side_fluid_temperature.get::<kelvin>()
                        )
                    );

                let wall_prandtl_number_part_correct: Ratio 
                    = fluid_material.try_get_prandtl_liquid(
                        part_correct_wall_temperature,
                        atmospheric_pressure
                    )?;

                pipe_prandtl_reynolds_gnielinksi_data.prandtl_wall = wall_prandtl_number_part_correct;

                nusselt_estimate_shell_side = shell_side_fluid_to_inner_tube_surf_nusselt_correlation.
                    estimate_based_on_prandtl_reynolds_and_wall_correction(
                        bulk_prandtl_number, 
                        wall_prandtl_number_part_correct,
                        reynolds_number_shell_side_abs_for_nusselt_estimate)?;

            } else {
                nusselt_estimate_shell_side = shell_side_fluid_to_inner_tube_surf_nusselt_correlation.
                    estimate_based_on_prandtl_and_reynolds_no_wall_correction(
                        bulk_prandtl_number, 
                        reynolds_number_shell_side_abs_for_nusselt_estimate)?;

            }



            // now we can get the heat transfer coeff, 

            let shell_side_h_to_fluid: HeatTransfer;

            let k_fluid_average: ThermalConductivity = 
                fluid_material.try_get_thermal_conductivity(
                    shell_side_fluid_temperature)?;

            shell_side_h_to_fluid = nusselt_estimate_shell_side * k_fluid_average / shell_side_fluid_hydraulic_diameter;

            //// this was for debugging
            //dbg!(
            //    &(nusselt_estimate_shell_side,
            //        nusselt_estimate_tube_side,
            //        shell_side_fluid_temperature.get::<uom::si::thermodynamic_temperature::degree_celsius>(),
            //        reynolds_number_shell_side_abs_for_nusselt_estimate,
            //        reynolds_number_single_tube_abs_for_nusselt_estimate)
            //);

            // 1/h_s 
            let one_over_hs = 1.0 as f64  / shell_side_h_to_fluid;



            // 1/U = 1/h_t d_o/d_i + d_o/(2 lambda_w) ln (d_o/d_i) 
            // + 1/h_s


            let one_over_u = 
                one_over_ht_times_do_by_di + 
                do_by_2_lambda_w_times_ln_do_by_di +
                one_over_hs;



            // overall heat transfer coeff shell side
            let u_based_on_shell_side_area: HeatTransfer = 
                1.0 as f64 / one_over_u;

            Ok(u_based_on_shell_side_area)
        }

    /// provides the tube bundle side heat transfer area 
    /// on the shell side
    pub fn tube_bundle_heat_transfer_area_shell_side(&self) -> Area {

        let n_t = self.number_of_tubes;
        let pipe_shell_clone: SolidColumn = 
            self.inner_pipe_shell_array_for_single_tube.
            clone().try_into().unwrap();

        let l = pipe_shell_clone.get_component_length();

        let d_o = self.tube_side_od;

        let shell_area: Area = 
            n_t as f64 * PI * d_o * l;

        return shell_area;


    }


}
