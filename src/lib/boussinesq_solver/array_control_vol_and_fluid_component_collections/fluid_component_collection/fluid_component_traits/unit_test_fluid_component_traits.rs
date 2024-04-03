
use uom::si::{dynamic_viscosity::poise, mass_density::kilogram_per_cubic_meter, mass_rate::kilogram_per_second, pressure::pascal, ratio::ratio};





#[test]
pub fn test_custom_fluid_component_pressure_loss() -> Result<(),
crate::thermal_hydraulics_error::ThermalHydraulicsLibError>{
    // this tests the calc pressure loss for fluid component 
    use uom::si::f64::*;
    use uom::ConstZero;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidCustomComponentCalcPressureLoss;
    use uom::si::length::{meter, inch, millimeter};
    use uom::si::mass_rate::kilogram_per_hour;
    use uom::si::angle::degree;
    use std::f64::consts::PI;

    // this is the test component
    pub struct TestComponent {
        pub loss_correlation: DimensionlessDarcyLossCorrelations
    }

    impl FluidCustomComponentCalcPressureLoss for TestComponent {
        fn get_custom_loss_correlations(&mut self) ->
            crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations {
                self.loss_correlation.clone()
        }

        fn get_custom_loss_correlations_immutable(&self) ->
            crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations {
                self.loss_correlation.clone()
        }

        fn set_custom_loss_correlations(
            &mut self,
            custom_loss_correlation: crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations) {
            self.loss_correlation = custom_loss_correlation;
        }

        fn get_custom_component_absolute_roughness(
            &mut self) -> Length {
            Length::ZERO
        }

        fn get_custom_component_absolute_roughness_immutable(
            &self) -> Length {
            Length::ZERO
        }
    }

    // let's have a test mass flowrate 

    let fluid_mass_flowrate = MassRate::new::<kilogram_per_second>(0.2);
    let hydraulic_diameter = Length::new::<inch>(1.0);
    let incline_angle = Angle::new::<degree>(90.0);
    let component_length = Length::new::<meter>(0.5);
    let absolute_roughness = Length::new::<millimeter>(0.001);
    let cross_sectional_area = hydraulic_diameter*hydraulic_diameter*PI/4.0_f64;
    let fluid_viscosity = DynamicViscosity::new::<poise>(0.01);
    let fluid_density = MassDensity::new::<kilogram_per_cubic_meter>(1000.0);

    let a = Ratio::new::<ratio>(18.0);
    let b = Ratio::new::<ratio>(93000_f64);
    let c: f64 = -1.35;

    let loss_correlation = DimensionlessDarcyLossCorrelations::new_simple_reynolds_power_component(
        a, b, c);
    
    // create object 
    let test_component_object = TestComponent{ loss_correlation };
    // calculate a pressure loss

    {
        // forward test 

        let pressure_loss: Pressure = TestComponent::fluid_custom_component_calc_pressure_loss(
            fluid_mass_flowrate, 
            cross_sectional_area, 
            hydraulic_diameter, 
            fluid_viscosity, 
            fluid_density, 
            test_component_object.loss_correlation).unwrap();

        // expected pressure loss is 1430 pascals
        approx::assert_relative_eq!(
            1430.0,
            pressure_loss.get::<pascal>(),
            max_relative=0.01);
    }

    {
        // reverse test 

        let pressure_loss: Pressure = TestComponent::fluid_custom_component_calc_pressure_loss(
            -fluid_mass_flowrate, 
            cross_sectional_area, 
            hydraulic_diameter, 
            fluid_viscosity, 
            fluid_density, 
            test_component_object.loss_correlation).unwrap();

        // expected pressure loss is -1430 pascals
        // (other direction)
        approx::assert_relative_eq!(
            -1430.0,
            pressure_loss.get::<pascal>(),
            max_relative=0.01);
    }


    Ok(())
}
