





/// Example 3,
/// 
/// suppose now we have a coriolis flowmeter
/// with a custom friction factor correlation
///
/// (f_darcy L/D + K) = 18 + 93000/Re^1.35
///
/// we shall use water to push flow through this coriolis flowmeter
///
/// also, the programming is rather tedious
/// because of lifetimes, but this is one example of how it can be done
#[test]
//#[ignore="debugging"]
pub fn coriolis_flowmeter_empirical_custom_component_example_3(){

    use std::f64::consts::PI;

    use uom::si::dynamic_viscosity::poise;
    use uom::si::f64::*;
    use uom::si::length::{meter, inch, millimeter};
    use uom::si::mass_density::kilogram_per_cubic_meter;
    use uom::si::mass_rate::kilogram_per_second;
    use uom::si::pressure::pascal;
    use uom::si::angle::degree;
    use uom::si::ratio::ratio;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidComponentTrait;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidCustomComponentCalcPressureChange;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::fluid_component_collection::fluid_component_traits::FluidCustomComponentCalcPressureLoss;
    use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::fluid_component_calculation::DimensionlessDarcyLossCorrelations;

    struct CoriolisFlowmeter {

        pressure_loss: Pressure,
        mass_flowrate: MassRate,
        internal_pressure: Pressure,
        hydraulic_diameter: Length,
        incline_angle: Angle,
        component_length: Length,
        fluid_density: MassDensity,
        fluid_viscosity: DynamicViscosity,
        absolute_roughness: Length,
        loss_correlation: DimensionlessDarcyLossCorrelations
    }

    impl FluidCustomComponentCalcPressureChange for CoriolisFlowmeter {

    }

    impl FluidCustomComponentCalcPressureLoss for CoriolisFlowmeter {

        fn get_custom_component_absolute_roughness(
            &mut self) -> Length {

            return self.absolute_roughness;
        }

        fn get_custom_component_absolute_roughness_immutable(
            &self) -> Length {

            return self.absolute_roughness;
        }

        fn get_custom_loss_correlations(&mut self) -> DimensionlessDarcyLossCorrelations {
            self.loss_correlation.clone()
        }

        fn get_custom_loss_correlations_immutable(&self) -> DimensionlessDarcyLossCorrelations {
            self.loss_correlation.clone()
        }

        fn set_custom_loss_correlations(
            &mut self,
            custom_loss_correlation: DimensionlessDarcyLossCorrelations) {
            self.loss_correlation = custom_loss_correlation;
        }




    }

    impl FluidComponentTrait for CoriolisFlowmeter {
        fn set_internal_pressure_source(
            &mut self,
            internal_pressure: Pressure) {

            self.internal_pressure = internal_pressure;
        }


        fn get_internal_pressure_source(
            &mut self) -> Pressure{
            return self.internal_pressure;
        }

        fn get_internal_pressure_source_immutable(
            &self) -> Pressure{
            return self.internal_pressure;
        }

        fn set_mass_flowrate(
            &mut self,
            mass_flowrate: MassRate){

            self.mass_flowrate = mass_flowrate;
        }

        fn set_pressure_loss(
            &mut self,
            pressure_loss: Pressure){
            self.pressure_loss = pressure_loss;
        }

        fn get_hydraulic_diameter(&mut self) -> Length{

            return self.hydraulic_diameter;

        }

        fn get_hydraulic_diameter_immutable(
            &self) -> Length{

            return self.hydraulic_diameter;

        }

        fn get_incline_angle(&mut self) -> Angle {

            return self.incline_angle;

        }

        fn get_incline_angle_immutable(&self) -> Angle {

            return self.incline_angle;

        }

        fn get_component_length(&mut self) -> Length {

            return self.component_length;
        }

        fn get_component_length_immutable(&self) -> Length {

            return self.component_length;
        }

        fn get_fluid_density_at_ref_temperature(&mut self) -> MassDensity {

            return self.fluid_density;

        }

        fn get_fluid_density_immutable_at_ref_temperature(&self) -> MassDensity {

            return self.fluid_density;

        }

        fn get_fluid_viscosity_at_ref_temperature(&mut self) -> DynamicViscosity {

            return self.fluid_viscosity;

        }

        fn get_fluid_viscosity_immutable_at_ref_temperature(&self) -> DynamicViscosity {

            return self.fluid_viscosity;

        }

        fn get_cross_sectional_area(&mut self) -> Area {

            return self.get_hydraulic_diameter()*
                self.get_hydraulic_diameter()*
                PI/4.0_f64;

        }

        fn get_cross_sectional_area_immutable(&self) -> Area {

            return self.get_hydraulic_diameter_immutable()*
                self.get_hydraulic_diameter_immutable()*
                PI/4.0_f64;

        }


        /// gets pressure loss given current state of
        /// the component 
        fn get_pressure_loss(&mut self) -> Pressure {

            let fluid_mass_flowrate = 
                self.mass_flowrate;

            let cross_sectional_area = 
                self.get_cross_sectional_area();

            let hydraulic_diameter = 
                self.get_hydraulic_diameter();

            let fluid_viscosity = 
                self.get_fluid_viscosity_at_ref_temperature();

            let fluid_density = 
                self.get_fluid_density_at_ref_temperature();


            let pressure_loss =
                CoriolisFlowmeter::fluid_custom_component_calc_pressure_loss(
                    fluid_mass_flowrate, 
                    cross_sectional_area, 
                    hydraulic_diameter, 
                    fluid_viscosity, 
                    fluid_density, 
                    self.loss_correlation
                    ).unwrap();

            self.pressure_loss = pressure_loss;

            return pressure_loss;

        }

        /// gets pressure loss given current state
        /// of the system except for mass flowrate
        /// with an immutable borrow of self
        fn get_pressure_loss_immutable(
            &self,
            mass_flowrate: MassRate) -> Pressure {

            let fluid_mass_flowrate = 
                mass_flowrate;

            let cross_sectional_area = 
                self.get_cross_sectional_area_immutable();

            let hydraulic_diameter = 
                self.get_hydraulic_diameter_immutable();

            let fluid_viscosity = 
                self.get_fluid_viscosity_immutable_at_ref_temperature();

            let fluid_density = 
                self.get_fluid_density_immutable_at_ref_temperature();


            // i need to make some immutable borrows here...

            let pressure_loss =
                CoriolisFlowmeter::fluid_custom_component_calc_pressure_loss(
                    fluid_mass_flowrate, 
                    cross_sectional_area, 
                    hydraulic_diameter, 
                    fluid_viscosity, 
                    fluid_density, 
                    self.loss_correlation, ).unwrap();


            return pressure_loss;

        }

        /// gets mass flowrate given current state
        /// of the pipe
        fn get_mass_flowrate(&mut self) -> MassRate {

            //i'll have to get the pressure change
            //
            // pressure_change = 
            // - pressure_change
            // + hydrostatic pressure change
            // + internal pressure source
            //

            // internal pressure source
            let internal_pressure_source = 
                self.get_internal_pressure_source();

            // hydrostatic pressure
            let incline_angle = 
                self.get_incline_angle();

            let hydrostatic_pressure_change =
                self.get_hydrostatic_pressure_change_at_ref_temperature();

            // pressure_loss term
            //
            //
            let pressure_loss = 
                self.get_pressure_loss();

            // now we get pressure change

            let pressure_change =
                - pressure_loss
                + hydrostatic_pressure_change
                + internal_pressure_source;


            let cross_sectional_area = 
                self.get_cross_sectional_area();

            let hydraulic_diameter = 
                self.get_hydraulic_diameter();

            let fluid_viscosity = 
                self.get_fluid_viscosity_at_ref_temperature();

            let fluid_density = 
                self.get_fluid_density_at_ref_temperature();


            let absolute_roughness = 
                self.get_custom_component_absolute_roughness();

            let source_pressure = 
                self.get_internal_pressure_source();

            let mass_flowrate =
                CoriolisFlowmeter::
                fluid_custom_component_calc_mass_flowrate_from_pressure_change(
                    pressure_change, 
                    cross_sectional_area, 
                    hydraulic_diameter, 
                    fluid_viscosity, 
                    fluid_density, 
                    absolute_roughness, 
                    incline_angle, 
                    source_pressure, 
                    self.loss_correlation).unwrap();

            self.mass_flowrate = mass_flowrate;

            return mass_flowrate;
        }

        /// gets mass flowrate given current state of the pipe
        /// except for pressure loss
        fn get_mass_flowrate_from_pressure_loss_immutable(
            &self,
            pressure_loss: Pressure) -> MassRate {

            //i'll have to get the pressure change
            //
            // pressure_change = 
            // - pressure_change
            // + hydrostatic pressure change
            // + internal pressure source
            //

            // internal pressure source
            let internal_pressure_source = 
                self.get_internal_pressure_source_immutable();

            // hydrostatic pressure

            let incline_angle = 
                self.get_incline_angle_immutable();


            let hydrostatic_pressure_change =
                self.get_hydrostatic_pressure_change_immutable_at_ref_temperature();


            // now we get pressure change

            let pressure_change =
                - pressure_loss
                + hydrostatic_pressure_change
                + internal_pressure_source;


            let cross_sectional_area = 
                self.get_cross_sectional_area_immutable();

            let hydraulic_diameter = 
                self.get_hydraulic_diameter_immutable();

            let fluid_viscosity = 
                self.get_fluid_viscosity_immutable_at_ref_temperature();

            let fluid_density = 
                self.get_fluid_density_immutable_at_ref_temperature();

            let absolute_roughness = 
                self.get_custom_component_absolute_roughness_immutable();

            let source_pressure = 
                self.get_internal_pressure_source_immutable();

            let mass_flowrate =
                CoriolisFlowmeter::
                fluid_custom_component_calc_mass_flowrate_from_pressure_change(
                    pressure_change, 
                    cross_sectional_area, 
                    hydraulic_diameter, 
                    fluid_viscosity, 
                    fluid_density, 
                    absolute_roughness, 
                    incline_angle, 
                    source_pressure, 
                    self.loss_correlation).unwrap();

            return mass_flowrate;
        }
    }

    impl CoriolisFlowmeter {

        /// Reynold's power correlation in the form 
        /// f_darcy = A + B Re^(C)
        ///
        /// The first in the tuple is A,
        /// the second is B, the third is C
        fn new(hydraulic_diameter: Length,
               incline_angle: Angle,
               component_length: Length,
               absolute_roughness: Length,
               a: Ratio,
               b: Ratio,
               c: f64 ) -> Self {

            let loss_correlation: DimensionlessDarcyLossCorrelations = 
                DimensionlessDarcyLossCorrelations::new_simple_reynolds_power_component(
                    a, b, c);

            // by default, i set pressure loss and mass flowrate to 0 
            // internal pressure also set to 0
            return Self { 
                pressure_loss: Pressure::new::<pascal>(0.0), 
                mass_flowrate: MassRate::new::<kilogram_per_second>(0.0), 
                internal_pressure: Pressure::new::<pascal>(0.0),
                hydraulic_diameter, 
                incline_angle, 
                component_length, 
                fluid_density: MassDensity::new::<kilogram_per_cubic_meter>(1000.0), 
                fluid_viscosity: DynamicViscosity::new::<poise>(0.01),
                absolute_roughness,
                loss_correlation,
            }


        }

    }

    // now we have defined our coriolis flowmeter with water, we can start!

    let hydraulic_diameter = 
        Length::new::<inch>(1.0);

    let incline_angle = 
        Angle::new::<degree>(90.0);

    let component_length = 
        Length::new::<meter>(0.5);

    let absolute_roughness = 
        Length::new::<millimeter>(0.001);


    let mut flowmeter_object = 
        CoriolisFlowmeter::new(
            hydraulic_diameter, 
            incline_angle, 
            component_length, 
            absolute_roughness, 
            Ratio::new::<ratio>(18.0),
            Ratio::new::<ratio>(93000_f64),
            -1.35);


    // set mass flowrate at 0.2 kg/s

    let mut mass_flowrate = MassRate::new::<kilogram_per_second>(0.2);

    flowmeter_object.set_mass_flowrate(mass_flowrate);

    let mut pressure_change = 
        flowmeter_object.get_pressure_change();


    // expected pressure loss is 1430 pascals
    // expected pressure change is -6335 pascals
    // becuase of elevation
    approx::assert_relative_eq!(
        -6335_f64,
        pressure_change.value,
        max_relative=0.01);

    // we'll now test the mass flowrate portion

    pressure_change = Pressure::new::<pascal>(-6335_f64);

    flowmeter_object.set_pressure_change(pressure_change);

    mass_flowrate = flowmeter_object.get_mass_flowrate();


    approx::assert_relative_eq!(
        0.2,
        mass_flowrate.value,
        max_relative=0.01);

    // of these functions and see if they work well
    //
    // the immutable versions of the methods take in &self rather
    // than &mut self, this enables safety in terms of parallelism
    // and may help with the use of peroxide iteration libraries
    // which are numerical root finders. These root finders 
    // cannot use mutable functions


    approx::assert_relative_eq!(
        mass_flowrate.value,
        flowmeter_object.
        get_mass_flowrate_from_pressure_change_immutable(pressure_change)
        .value,
        max_relative=0.01);

    // now we can get pressure loss in both direction
    // one should be the negative value of the other if
    // done correctly...

    flowmeter_object.set_mass_flowrate(
        MassRate::new::<kilogram_per_second>(0.2));

    let pressure_loss_positive_direction = 
        flowmeter_object.get_pressure_loss();

    flowmeter_object.set_mass_flowrate(
        MassRate::new::<kilogram_per_second>(-0.2));

    let pressure_loss_negative_direction = 
        flowmeter_object.get_pressure_loss();

    approx::assert_relative_eq!(
        pressure_loss_positive_direction.value,
        1430_f64,
        max_relative=0.01);

    approx::assert_relative_eq!(
        pressure_loss_positive_direction.value,
        -pressure_loss_negative_direction.value,
        max_relative=0.01);

}

