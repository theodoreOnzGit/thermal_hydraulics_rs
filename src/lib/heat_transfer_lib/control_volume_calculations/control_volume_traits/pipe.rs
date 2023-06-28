use std::f32::consts::PI;

use uom::{si::f64::*, num_traits::ToPrimitive};
use crate::heat_transfer_lib::nusselt_correlations
::pipe_correlations::gnielinski_correlation_interpolated_uniform_heat_flux_liquids_developing;

use super::{super::common_functions::*, ControlVolume};
use uom::si::power::watt;

/// contains traits specific to pipe flows
///
/// this means flows for one mass flowrate in,
/// one mass flowrate out
///
///
/// the pipe is constructed as follows.
///
/// It has a "front" and "back"
///
/// if flow is going forwards,
/// flow goes from "front" to "back"
///
/// if flow is going backwards
/// flow goes from "back" to "front"
trait PipeThermalComponent : ControlVolume{

    /// get pipe back temperature
    fn get_pipe_back_temperature(&self) -> ThermodynamicTemperature;

    /// get pipe front temperature
    fn get_pipe_front_temperature(&self) -> ThermodynamicTemperature;

    /// get mass flowrate through the pipe,
    /// positive means that it's forward flow
    /// negative means it's backward flow
    fn get_mass_flowrate_through_pipe(&self) -> MassRate;

    /// gets heat loss for pipe to surroundings
    fn get_pipe_heat_loss(&self) -> Power;

    /// get heat gain for pipe
    fn get_pipe_heat_gain(&self) -> Power;

    /// get pipe work done, normally zero
    fn get_pipe_work_done(&self) -> Power {
        return Power::new::<watt>(0.0);
    }


    /// get pipe enthalpy outflow
    fn get_pipe_enthalpy_outflow(&self) -> Power {

        // first we get front and back temperatures
        let front_temp = self.get_pipe_front_temperature();
        let back_temp = self.get_pipe_back_temperature();

        // i have a variable for outlet thermodynamic temp
        // no value yet, but settoutg type
        let outlet_temp: ThermodynamicTemperature;

        // then we check which one is outlet temperature
        if self.pipe_is_forward_flow() {
            outlet_temp = back_temp.clone();
        } else {
            outlet_temp = front_temp.clone();
        }

        // next we get the specific enthalpy
        let outlet_specific_enthalpy = 
            self.get_specific_fluid_enthalpy(outlet_temp);

        // we return the outlet enthalpy
        return outlet_specific_enthalpy * 
            self.get_mass_flowrate_through_pipe().abs();
    }

    /// get pipe enthalpy inflow
    fn get_pipe_enthalpy_inflow(&self) -> Power {

        // first we get front and back temperatures
        let front_temp = self.get_pipe_front_temperature();
        let back_temp = self.get_pipe_back_temperature();

        // i have a variable for inlet thermodynamic temp
        // no value yet, but setting type
        let inlet_temp: ThermodynamicTemperature;

        // then we check which one is inlet temperature
        if self.pipe_is_forward_flow() {
            inlet_temp = front_temp.clone();
        } else {
            inlet_temp = back_temp.clone();
        }

        // next we get the specific enthalpy
        let inlet_specific_enthalpy = 
            self.get_specific_fluid_enthalpy(inlet_temp);

        // we return the inlet enthalpy
        return inlet_specific_enthalpy * 
            self.get_mass_flowrate_through_pipe().abs();
    }


    /// check if forward flow is true for pipe
    fn pipe_is_forward_flow(&self) -> bool {

        let mass_flowrate = 
            self.get_mass_flowrate_through_pipe();

        if mass_flowrate.value < 0.0 {
            return false;
        } else {
            return true;
        }

    }

    /// convert temperature to enthalpy
    fn get_specific_fluid_enthalpy
        (&self,
         fluid_temperature: ThermodynamicTemperature) -> AvailableEnergy;

    

}

/// pipe heat loss associated functions
///
/// this trait contains functions useful for pipes
/// when thermal inertia is negligble for the pipe coverings
/// insulation and etc
///
/// This requires us to input thermal resistance parameters
/// for example thickness of pipe layers, 
/// ambient temperature and surrounding temperature
trait NoThermalResistancePipeHeatLoss{

    /// assumes the pipe is circular, and has
    /// two layers
    ///
    /// one layer which is the actual pipe, eg. copper
    /// and another layer is the thermal insulator
    /// eg. fiberglass or foam or something
    ///
    /// you are required to specify the inner and outer diameter of the pipe,
    /// insulation thickness,
    /// pipe bulk temperature
    /// entrance region length and etc
    /// 
    /// air heat convection coefficient is estimated at 20 W/(m^2 K)
    /// 
    /// you also need to provide a darcy friction factor as an argument
    fn calculate_heat_loss_rate_to_ambient_with_insulation(
        fluid_mass_flowrate: MassRate,
        inner_diameter: Length,
        fluid_bulk_temperature: ThermodynamicTemperature,
        wall_temperautre: ThermodynamicTemperature,
        darcy_friction_factor: f64,
        entrance_length_estimate: Length,
        pipe_length: Length,) -> 
        Result<Power,String> {

            let pi: f64 = PI.to_f64().unwrap();
            // let's first get the Re

            let inner_cross_sectional_area: Area
                = inner_diameter 
                * inner_diameter
                * pi
                / 4.0_f64;

            let fluid_viscosity = 
                Self::get_fluid_viscosity(fluid_bulk_temperature);

            let reynolds_number: Ratio = 
                fluid_mass_flowrate
                * inner_diameter
                / inner_cross_sectional_area
                / fluid_viscosity;

            // second, get prandtl number at wall and bulk

            let prandtl_wall: Ratio = 
                Self::get_fluid_prandtl_number(wall_temperautre);

            let prandtl_bulk: Ratio =
                Self::get_fluid_prandtl_number(fluid_bulk_temperature);

            // finally length to diameter ratio based on entrance length

            let length_to_diameter_ratio: Ratio =
                entrance_length_estimate
                / inner_diameter;

            // let's get our gnielinski nusselt number

            let nusselt_number_inner_wall =
                gnielinski_correlation_interpolated_uniform_heat_flux_liquids_developing(
                    reynolds_number.value, 
                    prandtl_bulk.value, 
                    prandtl_wall.value, 
                    darcy_friction_factor, 
                    length_to_diameter_ratio.value);

            // convert this into a heat transfer coefficient

            // Nu_D = hD/k
            // we take the D here to be inner diameter of the tube
            // secondly, we will need the thermal conductivity
            //

            let fluid_thermal_conductivity = 
                Self::get_fluid_thermal_conductivity(fluid_bulk_temperature);

            let inner_wall_heat_transfer_coefficient: HeatTransfer
                = fluid_thermal_conductivity 
                / inner_diameter
                * nusselt_number_inner_wall;
            
            // inner wall area pi D L

            let inner_wall_surface_area: Area 
                = inner_diameter
                * pipe_length
                * pi;

            // TODO:
            // 1) thermal resistance inputs for pipe wall 
            // 2) theraml resistance inputs for insulation
            // 3) thermal resistance for ambient air layer

            

            unimplemented!();
        }


    fn get_fluid_viscosity(fluid_temp: ThermodynamicTemperature) 
        -> DynamicViscosity;

    fn get_fluid_prandtl_number(fluid_temp: ThermodynamicTemperature,)
        -> Ratio;

    fn get_fluid_thermal_conductivity(fluid_temp: ThermodynamicTemperature)
        -> ThermalConductivity;
}

