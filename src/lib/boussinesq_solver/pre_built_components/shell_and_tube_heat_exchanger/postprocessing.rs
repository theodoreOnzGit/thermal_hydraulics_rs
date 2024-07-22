use super::SimpleShellAndTubeHeatExchanger;
use uom::si::{f64::*, pressure::atmosphere, ratio::ratio};
use crate::{boussinesq_solver::{array_control_vol_and_fluid_component_collections::{one_d_fluid_array_with_lateral_coupling::FluidArray, one_d_solid_array_with_lateral_coupling::SolidColumn}, boussinesq_thermophysical_properties::{LiquidMaterial, SolidMaterial}}, thermal_hydraulics_error::ThermalHydraulicsLibError};

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
            let reynolds_number_single_tube: Ratio = 
                single_tube_mass_flowrate/
                single_tube_flow_area
                *single_tube_hydraulic_diameter / viscosity;

            // next, bulk prandtl number 

            let bulk_prandtl_number: Ratio 
                = tube_fluid_material.try_get_prandtl_liquid(
                    tube_fluid_temperature,
                    atmospheric_pressure
                )?;

            let tube_side_nusselt_estimate: Ratio; 

            if correct_for_prandtl_wall_temperatures {

                let wall_prandtl_number: Ratio 
                    = tube_fluid_material.try_get_prandtl_liquid(
                        wall_temperature,
                        atmospheric_pressure
                    )?;

                tube_side_nusselt_estimate = 
                    self.tube_side_nusselt_correlation
                    .estimate_based_on_prandtl_reynolds_and_wall_correction(
                        bulk_prandtl_number, 
                        wall_prandtl_number, 
                        reynolds_number_single_tube)?;
            } else {
                tube_side_nusselt_estimate = 
                    self.tube_side_nusselt_correlation
                    .estimate_based_on_prandtl_reynolds_and_wall_correction(
                        bulk_prandtl_number, 
                        bulk_prandtl_number, 
                        reynolds_number_single_tube)?;
            }

            let k_fluid_average: ThermalConductivity = 
                tube_fluid_material.try_get_thermal_conductivity(
                    tube_fluid_temperature)?;

            let tube_side_h
                = tube_side_nusselt_estimate 
                * k_fluid_average / single_tube_hydraulic_diameter;
            
            let d_o_by_d_i: Ratio = 
                self.tube_side_od/self.tube_side_id;

            let ln_d_o_by_d_i: f64 = 
                d_o_by_d_i.get::<ratio>().ln();

            // to be continued...



            todo!();
    }

}
