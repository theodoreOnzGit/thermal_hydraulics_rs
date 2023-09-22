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
#[derive(Debug,Clone,PartialEq)]
pub struct HeaterTopBottomHead {

    inner_nodes: usize,

    pub twisted_tape_interior: HeatTransferEntity,

    pub steel_shell: HeatTransferEntity,

    pub therminol_array: HeatTransferEntity,

    pub ambient_temperature: ThermodynamicTemperature,

    pub heat_transfer_to_air: HeatTransfer,


}

impl HeaterTopBottomHead {

    /// synthesiszes a new heater top head
    ///
    /// K = 3.75 (this is in addition to pipe form losses)
    /// l = 0.0889m 
    /// d_h = 6.60e-3 m (hydraulic diameter)
    ///
    /// no insulation is assumed for this because 
    /// currently the ciet model has no insulation
    ///
    /// now, we can attach heat structures to either one of the control volumes 
    ///
    /// However, it should be attached to a heat structure 
    /// (perhaps a SteelSS304L single CV with appropriate conductance 
    /// 0.25 in^2 area, and callibrate length)
    ///
    /// It's only a rough guess
    ///
    /// We can callibrate until the outlet temperature at MX-10 is 102.45C 
    /// That was the steady state temperature not considering other 
    /// parasitic losses
    ///
    pub fn new_top_head(initial_temperature: ThermodynamicTemperature,
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
        // now twisted_tape 

        // the twisted tape width is assumed to be the twisted 
        // tape diameter in De Wet's dissertation
        let twisted_tape_width: Length = Length::new::<inch>(1.0);
        let twisted_tape_thickness = Length::new::<inch>(0.048);
        let twisted_tape_height = head_length;

        let twisted_tape = 
        SolidColumn::new_block(
            twisted_tape_height,
            twisted_tape_thickness,
            twisted_tape_width,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { inner_nodes: user_specified_inner_nodes,
            twisted_tape_interior: twisted_tape.into(),
            steel_shell: steel_shell_array.into(),
            therminol_array: therminol_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
        };
    }
    /// synthesiszes a new user callibrated heater top head
    ///
    /// K = 3.75 (this is in addition to pipe form losses)
    /// l = 0.0889m 
    /// d_h = 6.60e-3 m (hydraulic diameter)
    ///
    /// no insulation is assumed for this because 
    /// currently the ciet model has no insulation
    ///
    /// now, we can attach heat structures to either one of the control volumes 
    ///
    /// However, it should be attached to a heat structure 
    /// (perhaps a SteelSS304L single CV with appropriate conductance 
    /// 0.25 in^2 area, and callibrate length)
    ///
    /// It's only a rough guess
    ///
    /// We can callibrate until the outlet temperature at MX-10 is 102.45C 
    /// That was the steady state temperature not considering other 
    /// parasitic losses
    ///
    pub fn _new_user_callibrated_top_head(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature,
        h_to_air: HeatTransfer) -> Self {

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
        // now twisted_tape 

        // the twisted tape width is assumed to be the twisted 
        // tape diameter in De Wet's dissertation
        let twisted_tape_width: Length = Length::new::<inch>(1.0);
        let twisted_tape_thickness = Length::new::<inch>(0.048);
        let twisted_tape_height = head_length;

        let twisted_tape = 
        SolidColumn::new_block(
            twisted_tape_height,
            twisted_tape_thickness,
            twisted_tape_width,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { inner_nodes: user_specified_inner_nodes,
            twisted_tape_interior: twisted_tape.into(),
            steel_shell: steel_shell_array.into(),
            therminol_array: therminol_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
        };
    }

    /// synthesiszes a new heater bottom head
    ///
    /// K = 3.95 (this is in addition to pipe form losses)
    /// l = 0.19685m 
    /// d_h = 6.60e-3 m (hydraulic diameter)
    ///
    /// no insulation is assumed for this because 
    /// currently the ciet model has no insulation
    ///
    /// now, we can attach heat structures to either one of the control volumes 
    ///
    /// However, it should be attached to a heat structure 
    /// (perhaps a SteelSS304L single CV with appropriate conductance 
    /// 0.25 in^2 area, and callibrate length)
    ///
    /// It's only a rough guess
    ///
    pub fn new_bottom_head(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature) -> Self {

        let user_specified_inner_nodes: usize = 0;
        let flow_area = Area::new::<square_meter>(0.00105);
        let head_length = Length::new::<meter>(0.19685);
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);
        let pipe_form_loss = Ratio::new::<ratio>(3.95);

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
        let twisted_tape_width: Length = Length::new::<inch>(1.0);
        let twisted_tape_thickness = Length::new::<inch>(0.048);
        let twisted_tape_height = head_length;
        let twisted_tape = 
        SolidColumn::new_block(
            twisted_tape_height,
            twisted_tape_thickness,
            twisted_tape_width,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { inner_nodes: user_specified_inner_nodes,
            twisted_tape_interior: twisted_tape.into(),
            steel_shell: steel_shell_array.into(),
            therminol_array: therminol_array.into(),
            ambient_temperature,
            heat_transfer_to_air: h_to_air,
        };
    }
    /// synthesiszes a new heater bottom head allowing the 
    /// user to callibrate the heat transfer coefficient to air
    ///
    /// K = 3.95 (this is in addition to pipe form losses)
    /// l = 0.19685m 
    /// d_h = 6.60e-3 m (hydraulic diameter)
    ///
    /// no insulation is assumed for this because 
    /// currently the ciet model has no insulation
    ///
    /// now, we can attach heat structures to either one of the control volumes 
    ///
    /// However, it should be attached to a heat structure 
    /// (perhaps a SteelSS304L single CV with appropriate conductance 
    /// 0.25 in^2 area, and callibrate length)
    ///
    /// It's only a rough guess
    ///
    pub fn _new_user_callibrated_bottom_head(initial_temperature: ThermodynamicTemperature,
        ambient_temperature: ThermodynamicTemperature,
        h_to_air: HeatTransfer) -> Self {

        let user_specified_inner_nodes: usize = 0;
        let flow_area = Area::new::<square_meter>(0.00105);
        let head_length = Length::new::<meter>(0.19685);
        let atmospheric_pressure = Pressure::new::<atmosphere>(1.0);
        let pipe_form_loss = Ratio::new::<ratio>(3.95);

        // heater is inclined 90 degrees upwards, not that this is 
        // particularly important for this scenario

        let pipe_incline_angle = Angle::new::<uom::si::angle::degree>(90.0);

        // theoretically it's 6 W/(m^2 K) but then we'll have to manually 
        // input wall structures for additional heat loss
        //
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
        let twisted_tape_width: Length = Length::new::<inch>(1.0);
        let twisted_tape_thickness = Length::new::<inch>(0.048);
        let twisted_tape_height = head_length;
        let twisted_tape = 
        SolidColumn::new_block(
            twisted_tape_height,
            twisted_tape_thickness,
            twisted_tape_width,
            initial_temperature,
            atmospheric_pressure,
            SolidMaterial::SteelSS304L,
            user_specified_inner_nodes 
        );

        return Self { inner_nodes: user_specified_inner_nodes,
            twisted_tape_interior: twisted_tape.into(),
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
