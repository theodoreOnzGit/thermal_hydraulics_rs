use super::HeaterVersion2Bare;
use thermal_hydraulics_rs::prelude::alpha_nightly::*;
use uom::si::area::square_inch;

impl HeaterVersion2Bare {

    pub fn get_air_steel_shell_conductance(&mut self) 
        -> ThermalConductance {
        todo!()
    }

    pub fn get_therminol_node_steel_shell_conductance(&mut self) 
        -> ThermalConductance {
        todo!()
    }

    pub fn get_therminol_node_twisted_tape_conductance(
    &mut self) -> ThermalConductance {

        // the twisted tape itself acts as a thermal 
        // mass and effectively the twisted tape in the 
        // heater acts like porous media
        //
        // I need to find the nusselt correlation 
        // (otherwise use the Wakao Correlation) 
        //
        // and I also need the mass of the twisted tape overall 
        //
        // From De Wet's Dissertation, 
        // The heat transfer area for the twisted tape is 
        // 719 in^2 
        //
        // Also, the volume fraction of fluid in the 
        // original heater as compared to the 
        // fluid in the entire loop was approximately 3\% 
        // with heater v1.0, now, 
        //
        // There were two inserts tested in Lukas's 
        // conference paper 
        // First, a 51\% open perforated insert 
        // and a 23\% open one
        //
        // The 51\% open insert was used for heater v2.0 
        //
        // compared to the annular tube, it had a 157\%  
        // increase in residence time. Or for a constant 
        // mass or volumetric flow rate, a 157\% increase in volume 
        //
        // also, the volume fraction increased from 3\% of the 
        // loop to 8.1\% of the loop
        //
        // For heater v1.0, 
        // the flow volume is about 42.12 in^3 (the fluid volume height 
        // is 198 cm which includes the heater top and bottom heads)
        // whereas the flow volume in heater v2.0 is about 127 in^3
        //
        // Taking heater outer tube inner diameter of 1.5 in
        // and height of 78 in 
        // we get flow volume of 137.83 in^3
        //
        // Which means the twisted tape plus perforrated tube is about 
        // 10 in^3 of steel. We can use this to estimate the 
        // thermal inertia...
        //
        // I'm not too sure about the 157\% increase in residence 
        // time, 
        //
        // Okay, from Lukas's paper, the volume fraction of the 
        // fluid within the core compared to the loop increased from 
        // 3.3% to 8.1%, so this is about a 2.45 times as much as before 
        // assuming loop volume is reasonably large, this is close 
        // enough to the 3 times volume increase using the main fluid 
        // volume as a reference point 
        //
        // Hence, using the main fluid volume is right. And I think 
        // 10 in^3 of steel is reasonable for a thermal inertia 
        // measurement
        //
        // Thus fluid volume as modelled increases about 3 times
        // from v1.0 to v2.0, I wonder why we only have a 157\% increase 
        // in residence time. did the mass flowrate change?
        //
        // Apparently, the twisted tape height in fluid is 198cm 
        // which extends beyond the heated sections
        // 
        // So, heat transfer area is 719 in^2 including heater heads 
        // and the twisted tube is about 78 in or 198 cm long, which 
        // includes both heater heads, 
        //
        // we can scale heat transfer area accordingly using heated 
        // length of about 163 cm
        // 
        // for nusselt number, it seems best to use the Wakao 
        // Correlation as that is suitable for pebble beds anyhow
        // it's a best estimate, not need to be perfect for now 
        //
        // Can't really do much until someone does a separate 
        // effects test (SET)
        // 

        // find suitable heat transfer area
        let heated_length = Length::new::<inch>(66.0);
        let heated_length_plus_heads = Length::new::<inch>(78.0);

        let heat_transfer_area_heated_length_plus_heads: Area = 
        Area::new::<square_inch>(719.0);

        let heat_transfer_area_heated_length_only: Area
        = heated_length/ heated_length_plus_heads * 
        heat_transfer_area_heated_length_plus_heads;

        // next, need the nusselt number based on Wakao Correlation 
        let mass_flowrate = self.get_mass_flowrate();
        let flow_area: Area = Area::new::<square_inch>(1.63);
        let viscosity = self.get_fluid_viscosity();
        let hydraulic_diameter = Length::new::<inch>(0.5776);

        // need to convert hydraulic diameter to an equivalent 
        // spherical diameter
        //
        // but for now, I'm going to use Re and Nu using hydraulic diameter 
        // and live with it for the time being
        //
        let reynolds: Ratio = 
        mass_flowrate/flow_area*hydraulic_diameter / viscosity;
        
        // need to get prandtl number of fluid 
        // so I need fluid temperature 
        let mut therminol_fluid_array: FluidArray = 
        self.therminol_array.clone().try_into().unwrap();

        let fluid_average_temperature: ThermodynamicTemperature 
        = therminol_fluid_array.get_bulk_temperature().unwrap();

        let fluid_average_pressure: Pressure 
        = therminol_fluid_array.pressure_control_volume;

        let fluid_material: LiquidMaterial = 
        therminol_fluid_array.material_control_volume.try_into().unwrap();

        let fluid_prandtl: Ratio = fluid_material.try_get_prandtl_liquid
            (fluid_average_temperature, fluid_average_pressure)
            .unwrap();

        

        todo!()
    }
}


