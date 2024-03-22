// marked for deprecation, no longer used
#![warn(missing_docs)]
extern crate uom;
use uom::si::available_energy::joule_per_kilogram;
use uom::si::f64::*;
use uom::si::*;
use uom::si::heat_capacity::joule_per_kelvin;
use crate::heat_transfer_lib::control_volume_calculations::ExplicitCalculations::*;
use crate::ControlVolumeCalculations::
FluidEntity_StructsAndTraits::*;


// here are imports for units
use uom::si::time::second;
use uom::si::volume::cubic_meter;
use uom::si::thermodynamic_temperature::kelvin;

extern crate fluid_mechanics_rust;
use fluid_mechanics_rust::therminol_component::
FluidProperties;
use fluid_mechanics_rust::therminol_component::
dowtherm_a_properties;


/// Contains functions which return the
/// viscosity, density, enthalpy, specific heat capacity
/// and thermal conductivity for therminol VP 1 or
/// dowtherm A in the range 20-180C
///
/// The dowtherm A correlations are used
/// 
pub trait TherminolFluidProperties {

    /// returns dowtherm A density given a temperature
    fn density(fluid_temp: ThermodynamicTemperature) 
        -> MassDensity {
        return dowtherm_a_properties::getDowthermADensity(fluid_temp);
    }

    /// returns dowtherm A dynamic viscosity given
    /// a temperature
    fn viscosity(
        fluid_temp: ThermodynamicTemperature) -> DynamicViscosity{
        return dowtherm_a_properties::getDowthermAViscosity(fluid_temp);
    }

    ///returns dowtherm A specific
    ///enthalpy  given a temperature
    ///
    /// 0 J/kg specific enthalpy is assumed at 20C
    /// and everything is calculated from there
    fn enthalpy(fluid_temp: ThermodynamicTemperature) -> AvailableEnergy{
        return dowtherm_a_properties::getDowthermAEnthalpy(fluid_temp);
    }

    /// returns dowtherm A specific heat capacity
    ///
    fn specific_heat_capacity(
        fluid_temp: ThermodynamicTemperature) -> SpecificHeatCapacity{
        return dowtherm_a_properties::
            getDowthermAConstantPressureSpecificHeatCapacity(
            fluid_temp);
    }

    /// returns dowtherm A thermal conductivity
    fn thermal_conductivity(
        fluid_temp: ThermodynamicTemperature) -> ThermalConductivity{
        return dowtherm_a_properties::
            getDowthermAThermalConductivity(fluid_temp);
    }

    ///returns dowtherm A temperature
    ///  given a specific enthalpy
    ///
    /// 0 J/kg specific enthalpy is assumed at 20C
    /// and everything is calculated from there
    fn get_temperature_from_enthalpy(
        fluid_enthalpy: AvailableEnergy) -> ThermodynamicTemperature{
        return dowtherm_a_properties::
            get_temperature_from_enthalpy(fluid_enthalpy);
    }

}

