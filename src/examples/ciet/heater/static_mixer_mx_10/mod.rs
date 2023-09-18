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

    pub insulation: HeatTransferEntity,

    pub steel_shell: HeatTransferEntity,

    pub therminol_array: HeatTransferEntity,

    pub ambient_temperature: ThermodynamicTemperature,

    pub heat_transfer_to_air: HeatTransfer,


}

impl StaticMixerMX10 {

    pub fn new_static_mixer(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature) -> Self {

        let user_specified_inner_nodes: usize = 0;
        let flow_area = Area::new::<square_meter>(0.00105);
        let head_length = Length::new::<meter>(0.0889);
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);
        let pipe_form_loss = Ratio::new::<ratio>(3.75);

        // heater is inclined 90 degrees upwards, not that this is 
        // particularly important for this scenario

        let pipe_incline_angle = Angle::new::<uom::si::angle::degree>(90.0);

        // theoretically it's 6 W/(m^2 K) but then we'll have to manually 
        // input wall structures for additional heat loss
        //
        let h_to_air: HeatTransfer = 
        HeatTransfer::new::<watt_per_square_meter_kelvin>(6.0);
        let id = Length::new::<meter>(0.0381);
        let od = Length::new::<meter>(0.04);


        // inner therminol array
        let therminol_array: FluidArray = 
        FluidArray::new_odd_shaped_pipe(
            head_length,
            flow_area,
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
            head_length,
            id,
            od,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );
        // now twisted_tape (TBC)
        let insulation = 
        SolidColumn::new_cylindrical_shell(
            head_length,
            id,
            od,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { inner_nodes: user_specified_inner_nodes,
            insulation: insulation.into(),
            steel_shell: steel_shell_array.into(),
            therminol_array: therminol_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
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
