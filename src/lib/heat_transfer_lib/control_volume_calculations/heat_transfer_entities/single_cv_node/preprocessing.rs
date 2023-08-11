use crate::heat_transfer_lib::thermophysical_properties::{specific_enthalpy::temperature_from_specific_enthalpy, thermal_diffusivity::thermal_diffusivity, specific_heat_capacity::specific_heat_capacity, Material};

use super::SingleCVNode;
use uom::si::{f64::*, length::meter, power::watt, time::second};


impl SingleCVNode {

    /// this is a function to determine the relevant time scales 
    /// this one is based on conduction, which calculates timescales 
    /// based on mesh fourier number
    /// for stable conduction
    ///
    /// for an uneven volume, it will just take the shortest of 
    /// these lengthscales to determine a proper time scale
    #[inline]
    pub fn calculate_conduction_timestep(&self) -> Result<Time,String>{

        // first let us get relevant length scales 
        // for shortest time step, we will get the shortest length 
        // scale

        let lengthscale_stability_vec = 
        self.mesh_stability_lengthscale_vector.clone();


        // initiate a simple loop to find the shortest length scale
        let mut min_lengthscale = Length::new::<meter>(0.0);

        for length in lengthscale_stability_vec.iter() {
            if *length > min_lengthscale {
                min_lengthscale = *length;
            }
        }

        // now we can determine the conduction time scale regardless 
        // of whether it's solid or liquid (no gas implementations yet, 
        // this is alpha code, probably needs a few rewrites)
        //

        let max_mesh_fourier_number: f64 = 0.25;

        let control_vol_material = self.material_control_volume.clone();
        let control_vol_pressure = self.pressure_control_volume.clone();
        let cv_temperature = temperature_from_specific_enthalpy(
            self.material_control_volume, 
            self.current_timestep_control_volume_specific_enthalpy, 
            self.pressure_control_volume)?;


        let thermal_diffusivity_coeff: DiffusionCoefficient = 
        thermal_diffusivity(control_vol_material, 
            cv_temperature, 
            control_vol_pressure)?;

        // now let's obtain the minimum timescale for conduction

        let min_conduction_timescale = max_mesh_fourier_number * 
        min_lengthscale *
        min_lengthscale / 
        thermal_diffusivity_coeff;


        return Ok(min_conduction_timescale);
    }

    /// appends timestep constrained to fourier number stability
    #[inline]
    pub fn append_conduction_mesh_stability_timestep(&mut self, 
        lengthscale: Length) -> Result<Time, String> {

        // now we can determine the conduction time scale regardless 
        // of whether it's solid or liquid (no gas implementations yet, 
        // this is alpha code, probably needs a few rewrites)
        //

        let max_mesh_fourier_number: f64 = 0.25;

        let control_vol_material = self.material_control_volume.clone();
        let control_vol_pressure = self.pressure_control_volume.clone();
        let cv_temperature = temperature_from_specific_enthalpy(
            self.material_control_volume, 
            self.current_timestep_control_volume_specific_enthalpy, 
            self.pressure_control_volume)?;


        let thermal_diffusivity_coeff: DiffusionCoefficient = 
        thermal_diffusivity(control_vol_material, 
            cv_temperature, 
            control_vol_pressure)?;

        // now let's obtain the minimum timescale for conduction

        let min_conduction_timescale = max_mesh_fourier_number * 
        lengthscale *
        lengthscale / 
        thermal_diffusivity_coeff;


        return Ok(min_conduction_timescale);
    }
    

    /// calculates a time scale based on maximum temperature change 
    /// within one time step, 
    #[inline]
    fn get_max_timestep_based_on_max_temperature_change(
        &self,
        max_temperature_change: TemperatureInterval) 
    -> Result<Time,String> {

        // let's get the net power change first 

        // first thing is to sum up all enthalpy rate of change 
        let mut total_enthalpy_rate_change = 
        Power::new::<watt>(0.0);

        for enthalpy_chg_rate in 
            self.rate_enthalpy_change_vector.clone().iter() {

                total_enthalpy_rate_change += *enthalpy_chg_rate;
            }

        // we have Q = m c_p Delta T
        // Q = power * time 
        //
        // power * time = m c_p  Delta T
        //
        // time = m c_p Delta T / power 
        // for small temperature changes, assume cp constant

        let cv_mass_clone = self.mass_control_volume.clone();
        let cv_material = self.material_control_volume.clone();
        let cv_temperature = self.get_temperature()?;
        let cv_pressure = self.pressure_control_volume.clone();

        let cv_heat_capacity = specific_heat_capacity(
            cv_material,
            cv_temperature,
            cv_pressure,
        )?;

        // now we can calculate time step

        let time_step_based_on_max_temperature_change = 
        cv_mass_clone * cv_heat_capacity * max_temperature_change / 
        total_enthalpy_rate_change;

        // now return it to the environment 

        return Ok(time_step_based_on_max_temperature_change);

    }

    /// compiles a list of time steps based on various criteria, 
    ///
    ///
    #[inline]
    pub fn get_max_timestep(&mut self, 
        max_temperature_change: TemperatureInterval) -> Result<Time, String> {

        // let's calculate conduction time step first, based on 
        // thermal diffusivity and fourier number
        let conduction_timestep = self.calculate_conduction_timestep()?;

        // next I want to push it to the vector 

        self.max_timestep_vector.push(conduction_timestep);

        // let's match the material to solid or liquid 

        let cv_material = self.material_control_volume.clone();

        match cv_material {
            Material::Solid(_) => (),
            Material::Liquid(_) => todo!("need to calculate liquid timestep"),
        }

        // last but not least, check if temperature changes exceed a 
        // certain mount, recommend 0.5 C or 1 C for max temperature 
        // change in one time step
        //
        // inspired by GeN-Foam and other papers.
        
        let max_temperature_change_timestep: Time = 
        self.get_max_timestep_based_on_max_temperature_change(
            max_temperature_change)?;

        self.max_timestep_vector.push(max_temperature_change_timestep);

        // we still need the timestep changes from the control volumes 
        // linked to this cv.


        // initiate a simple loop to find the shortest time scale
        // that is "safe" out of all time steps
        let mut min_timescale = Time::new::<second>(75_f64);

        for time in self.max_timestep_vector.iter() {
            if *time < min_timescale {
                min_timescale = *time;
            }
        }
        return Ok(min_timescale);

    }
}
