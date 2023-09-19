use thermal_hydraulics_rs::prelude::alpha_nightly::*;

use uom::si::area::square_meter;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::ratio::ratio;
use uom::si::pressure::atmosphere;
/// represents heater version 2 without insulation 
/// This is because during 2018-ish, the heater insulation 
/// got burnt off and a lot of frequency response tests were done 
/// with insulation removed
///
/// Heater version 2 bare has no insulation
/// but it has a twisted tape interior
///
///
/// note that it only contains the heated section, not the top nor 
/// bottom heads
pub struct StaticMixerMX10 {

    inner_nodes: usize,

    pub insulation_array: HeatTransferEntity,

    pub steel_shell: HeatTransferEntity,

    pub therminol_array: HeatTransferEntity,

    pub ambient_temperature: ThermodynamicTemperature,

    pub heat_transfer_to_air: HeatTransfer,

    tube_inner_diameter: Length,

    tube_outer_diameter: Length, 

    insulation_inner_diameter: Length,

    insulation_outer_diameter: Length,

    flow_area: Area,

    darcy_loss_correlation: DimensionlessDarcyLossCorrelations,
}

impl StaticMixerMX10 {


    /// constructs the static mixer using the RELAP/SAM model 
    /// as a basis 
    ///
    /// length = 0.33 m 
    /// d_h = 2.79e-2
    /// Insulation thickness: 5.08 cm
    /// (fiberglass)
    /// number of nodes (including two ends): 2
    ///
    /// Nusselt Number Correlation: same as heater (assumed)
    /// because there is quite a lot of mixing going on
    /// within the mixer
    ///
    /// Reynolds Number Correlation: unimplemented!()
    ///
    ///
    /// Unheated Structure Thermal Inertia: ignored
    pub fn new_static_mixer(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature) -> Self {

        let user_specified_inner_nodes: usize = 0;
        let flow_area = Area::new::<square_meter>(6.11e-4);
        let component_length = Length::new::<meter>(0.33);
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);
        let hydraulic_diameter = Length::new::<meter>(2.79e-2);

        // form losses not implemented yet
        let pipe_form_loss = Ratio::new::<ratio>(3.75);

        // heater is inclined 90 degrees upwards, not that this is 
        // particularly important for this scenario

        let pipe_incline_angle = Angle::new::<uom::si::angle::degree>(90.0);

        // theoretically it's 6 W/(m^2 K) but then we'll have to manually 
        // input wall structures for additional heat loss
        //
        let h_to_air: HeatTransfer = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(6.0);

        let fiberglass_thickness = Length::new::<meter>(0.0508);

        let steel_id = Length::new::<meter>(0.0381);
        let steel_od = Length::new::<meter>(0.04);
        let fiberglass_id = steel_od;
        let fiberglass_od = fiberglass_id + 
        fiberglass_thickness + fiberglass_thickness;


        // inner therminol array
        let therminol_array: FluidArray = 
        FluidArray::new_cylinder(
            component_length,
            hydraulic_diameter,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            LiquidMaterial::TherminolVP1,
            pipe_form_loss,
            user_specified_inner_nodes,
            pipe_incline_angle
        );
        // now the outer steel array
        let steel_shell_array = 
        SolidColumn::new_cylindrical_shell(
            component_length,
            steel_id,
            steel_od,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );
        // insulation
        let insulation = 
        SolidColumn::new_cylindrical_shell(
            component_length,
            fiberglass_id,
            fiberglass_od,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::Fiberglass,
            user_specified_inner_nodes 
        );

        // f + L/D K = 21 + 4000/Re
        let darcy_loss_correlation = 
        DimensionlessDarcyLossCorrelations::
            new_simple_reynolds_power_component(
                Ratio::new::<ratio>(21.0),
                Ratio::new::<ratio>(4000.0),
                -1.0
            );

        return Self { inner_nodes: user_specified_inner_nodes,
            insulation_array: insulation.into(),
            steel_shell: steel_shell_array.into(),
            therminol_array: therminol_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
            tube_inner_diameter: steel_id,
            tube_outer_diameter: steel_od,
            insulation_inner_diameter: fiberglass_id,
            insulation_outer_diameter: fiberglass_od,
            flow_area,
            darcy_loss_correlation,
        };
    }

}




pub mod preprocessing;
pub use preprocessing::*;

pub mod fluid_entity;
pub use fluid_entity::*;


pub mod calculation;
pub use calculation::*;

pub mod postprocessing;
pub use postprocessing::*;
