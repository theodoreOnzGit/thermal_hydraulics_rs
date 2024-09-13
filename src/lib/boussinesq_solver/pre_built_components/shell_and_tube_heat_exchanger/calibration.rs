use std::f64::consts::PI;

use uom::si::{f64::*, ratio::ratio};

use crate::{prelude::beta_testing::{SolidColumn, SolidMaterial}, thermal_hydraulics_error::ThermalHydraulicsLibError};

use super::SimpleShellAndTubeHeatExchanger;

impl SimpleShellAndTubeHeatExchanger {

    /// assuming sthe insulation is cylindrical,
    /// get the thermal resistance
    ///
    /// ln (d_o/d_i) * 1/(2 pi L lambda_insulation)
    pub fn get_insulation_cylindrical_thermal_resistance(&self) -> ThermalResistance {

        let insulation_id = self.shell_side_od;
        let insulation_od = insulation_id + 2.0*self.insulation_thickness;

        let l = self.get_effective_length();

        let lambda_insulation: ThermalConductivity;

        let mut insulation_array_clone: SolidColumn = self.insulation_array
            .clone()
            .try_into()
            .unwrap();

        let insulation_array_temp: ThermodynamicTemperature = 
            insulation_array_clone
            .try_get_bulk_temperature()
            .unwrap();

        let insulation_array_material: SolidMaterial = 
            insulation_array_clone.material_control_volume
            .clone()
            .try_into()
            .unwrap();

        lambda_insulation = 
            insulation_array_material
            .try_get_thermal_conductivity(
                insulation_array_temp
            ).unwrap();


        let log_numerator_term: f64 = 
            (insulation_od/insulation_id)
            .get::<ratio>()
            .ln();

        let denominator_term: ThermalConductance = 
            2.0 * PI * lambda_insulation * l;


        return log_numerator_term/denominator_term;

    }


    /// assuming sthe outer shell is cylindrical,
    /// get the thermal resistance
    ///
    /// ln (d_o/d_i) * 1/(2 pi L lambda_insulation)
    pub fn get_outer_shell_cylindrical_thermal_resistance(&self) -> 
        Result<ThermalResistance, ThermalHydraulicsLibError> {


            if !self.heat_exchanger_has_insulation {

                // thermal resistance is not implemented if we don't have 
                // insulation,
                // probably want to have a nicer error for this in future
                todo!();
            }

            let outer_shell_id = self.shell_side_id;
            let outer_shell_od = self.shell_side_od;

            let l = self.get_effective_length();

            let lambda_outer_shell: ThermalConductivity;

            let mut outer_shell_array_clone: SolidColumn = self.outer_shell
                .clone()
                .try_into()
                .unwrap();

            let outer_shell_array_temp: ThermodynamicTemperature = 
                outer_shell_array_clone
                .try_get_bulk_temperature()
                .unwrap();

            let outer_shell_array_material: SolidMaterial = 
                outer_shell_array_clone.material_control_volume
                .clone()
                .try_into()
                .unwrap();

            lambda_outer_shell = 
                outer_shell_array_material
                .try_get_thermal_conductivity(
                    outer_shell_array_temp
                ).unwrap();


            let log_numerator_term: f64 = 
                (outer_shell_od/outer_shell_id)
                .get::<ratio>()
                .ln();

            let denominator_term: ThermalConductance = 
                2.0 * PI * lambda_outer_shell * l;


            return Ok(log_numerator_term/denominator_term);

    }

    /// assuming sthe inner tubes are cylindrical parallel tubes,
    /// get the thermal resistance
    ///
    /// ln (d_o/d_i) * 1/(2 pi L lambda_insulation N_t)
    pub fn get_inner_tubes_cylindrical_thermal_resistance(&self) -> 
        Result<ThermalResistance, ThermalHydraulicsLibError> {



            let inner_tube_id = self.tube_side_id;
            let inner_tube_od = self.tube_side_od;

            let l = self.get_effective_length();

            let lambda_inner_tubes: ThermalConductivity;

            let mut inner_single_tube_array_clone: SolidColumn = self.inner_pipe_shell_array_for_single_tube
                .clone()
                .try_into()
                .unwrap();

            let inner_tube_array_temp: ThermodynamicTemperature = 
                inner_single_tube_array_clone
                .try_get_bulk_temperature()
                .unwrap();

            let inner_tube_array_material: SolidMaterial = 
                inner_single_tube_array_clone.material_control_volume
                .clone()
                .try_into()
                .unwrap();

            lambda_inner_tubes = 
                inner_tube_array_material
                .try_get_thermal_conductivity(
                    inner_tube_array_temp
                ).unwrap();


            let log_numerator_term: f64 = 
                (inner_tube_od/inner_tube_id)
                .get::<ratio>()
                .ln();

            let denominator_term: ThermalConductance = 
                2.0 * PI * lambda_inner_tubes * l * self.number_of_tubes as f64;


            return Ok(log_numerator_term/denominator_term);

    }
}
