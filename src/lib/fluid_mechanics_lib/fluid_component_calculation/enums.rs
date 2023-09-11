use crate::fluid_mechanics_lib::churchill_friction_factor;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use roots::*;
use uom::{si::{f64::*, ratio::ratio}, num_traits::Zero, typenum::P2};

/// contains form loss or minor loss correlations for use 
///
/// This will return a friction factor if one wishes it
#[derive(Clone, PartialEq, Copy, Debug)]
pub enum DimensionlessDarcyLossCorrelations {
    /// standard pipe loss, must input 
    /// roughness ratio 
    ///
    /// and also a K ratio for generic form losses
    Pipe(Ratio, Ratio, Ratio),
    /// Reynold's power correlation in the form 
    /// f_darcy = A + B Re^(C)
    ///
    /// The first in the tuple is A,
    /// the second is B, the third is C
    SimpleReynoldsPower(Ratio, Ratio, f64),

    /// Ergun Equation 
    /// Ergun, S., & Orning, A. A. (1949). Fluid flow through 
    /// randomly packed columns and fluidized beds. Industrial 
    /// & Engineering Chemistry, 41(6), 1179-1184.
    Ergun
}

impl Default for DimensionlessDarcyLossCorrelations {
    /// the default is just K = 1
    fn default() -> Self {
        return Self::new_simple_reynolds_power_component (
            Ratio::new::<ratio>(1.0),
            Ratio::new::<ratio>(0.0),
            0.0,
        );
    }
}

impl DimensionlessDarcyLossCorrelations {

    /// creates a new pipe object
    pub fn new_pipe(
        pipe_length: Length,
        surface_roughness: Length,
        hydraulic_diameter: Length,
        form_loss: Ratio
    ) -> Self {
        let length_to_diameter_ratio: Ratio = 
        pipe_length/hydraulic_diameter;

        let roughness_ratio: Ratio = 
        surface_roughness/hydraulic_diameter;

        return Self::Pipe(
            roughness_ratio,
            length_to_diameter_ratio,
            form_loss);
    }

    /// creates a new simple reynolds power correlation object 
    /// in the form
    /// Reynold's power correlation in the form 
    /// f_darcy = A + B Re^(C)
    pub fn new_simple_reynolds_power_component(
        a: Ratio,
        b: Ratio,
        c: f64) -> Self {

        return Self::SimpleReynoldsPower(
            a, b, c);
    }

    /// gets the darcy friction factor based on reynolds number and 
    /// other fluid component properties
    ///
    /// the convention is to disregard directionality,
    /// so reverse flow will also return a positive friction_factor
    /// value
    #[inline]
    pub fn darcy_friction_factor_fldk(&self, reynolds_input: Ratio) -> 
    Result<Ratio, ThermalHydraulicsLibError> {

        // check for reverse flow
        let mut reverse_flow = false;
        if reynolds_input.value < 0.0 {
            reverse_flow = true;
        }

        let reynolds: Ratio;
        if reverse_flow {
            reynolds = reynolds_input * -1.0;
        } else {
            reynolds = reynolds_input;
        }

        // check for zero flow 
        
        if reynolds_input.is_zero() {
            // return zero friction factor
            // to ensure that pressure losses are zero

            return Ok(Ratio::new::<ratio>(0.0))
        }


        let darcy_value: f64 = match self {
            // uses churchill_friction_factor here
            DimensionlessDarcyLossCorrelations::Pipe(roughness_ratio,
                length_to_diameter,
                form_loss) => {

                    // total friction factor is 
                    // = (f L/D + K)
                    //
                    let total_friction_factor = 
                    churchill_friction_factor::darcy(
                        reynolds.get::<ratio>(),
                        roughness_ratio.get::<ratio>())
                        * length_to_diameter.get::<ratio>()
                    + form_loss.get::<ratio>();

                    total_friction_factor
            },
            // f L/D + K = A + B Re^(C)
            DimensionlessDarcyLossCorrelations::SimpleReynoldsPower(a, b, c) => {
                let friction_factor  = a.get::<ratio>() + b.get::<ratio>()
                * reynolds.get::<ratio>().powf(*c);
                
                friction_factor

            },
            DimensionlessDarcyLossCorrelations::Ergun => {
                todo!()
            },
        };

        Ok(Ratio::new::<ratio>(darcy_value))
    }


    /// obtains bejan number given a reynolds number 
    /// this time, we consider directionality, so if reynolds number is 
    /// negative, pressure loss is also negative
    #[inline]
    pub fn get_bejan_number_from_reynolds(&self, reynolds_input: Ratio,)
    -> Result<Ratio, ThermalHydraulicsLibError>{


        // this is the fldk term
        // it will take care of Re = 0 but not directionality
        let total_losses_coeff = self.darcy_friction_factor_fldk(reynolds_input)?;

        // bejan number is 0.5 * fldk * Re^2
        //
        // Re^2 is always positive, and fldk should be greater or 
        // equal to positive 

        let mut bejan_number = 0.5 * total_losses_coeff * 
            reynolds_input * reynolds_input;

        // be mindful of reverse flow
        let mut reverse_flow = false;

        if reynolds_input.value < 0.0 {
            reverse_flow = true;
        }

        if reverse_flow {
            bejan_number = -1.0 * bejan_number;
            return Ok(bejan_number);
        }

        return Ok(bejan_number);

    }

    /// obtains a reynolds number from a given bejan number 
    ///
    /// needs testing
    #[inline]
    pub fn get_reynolds_number_from_bejan(&self,
        bejan_input: Ratio) -> Result<Ratio,ThermalHydraulicsLibError>{

        // we have to make a pressure drop root 

        // first we need limits for maximum and minimum reynolds 
        // number 
        // by default, 1e12 should be enough 

        let reynolds_max_limit_abs = Ratio::new::<ratio>(1.0e12);
        let upper_limit: f64 = reynolds_max_limit_abs.get::<ratio>();
        let lower_limit: f64 = -upper_limit;

        let pressure_drop_root = |reynolds: f64| -> f64 {
            // i'm solving for
            // Be - 0.5*fLDK*Re^2 = 0 
            // the fLDK term can be calculated using

            let lhs_bejan = bejan_input;
            let rhs_bejan: Ratio = self.get_bejan_number_from_reynolds(reynolds.into())
            .unwrap();

            return lhs_bejan.get::<ratio>() - rhs_bejan.get::<ratio>();
        };



        let mut convergency = SimpleConvergency { eps:1e-8f64, max_iter:30 };

        let reynolds_number_result
        = find_root_brent(upper_limit,
            lower_limit,
            pressure_drop_root,
            &mut convergency
        );
        
        let reynolds_number:f64 = reynolds_number_result.unwrap();
        

        return Ok(Ratio::new::<ratio>(reynolds_number));
    }

    /// pressure drop from Re
    /// characteristic lengthscales and fluid properties
    #[inline]
    pub fn get_pressure_loss_from_reynolds(&self,
        reynolds_input: Ratio,
        hydraulic_diameter: Length,
        fluid_density: MassDensity,
        fluid_viscosity:DynamicViscosity) -> 
    Result<Pressure,ThermalHydraulicsLibError>{

        if fluid_viscosity.value <= 0.0 {
            panic!("fluid Viscosity <= 0.0, nonphysical");
        }

        if hydraulic_diameter.value <= 0.0 {
            panic!("hydraulic Diameter <= 0.0, nonphysical");
        }

        if fluid_density.value <= 0.0 {
            panic!("fluidDensity <= 0.0, nonphysical");
        }

        let bejan_number = self.get_bejan_number_from_reynolds(
            reynolds_input)?;

        let fluid_pressure = fluid_viscosity.powi(P2::new())*
        bejan_number/
        hydraulic_diameter.powi(P2::new())/
        fluid_density;

        return Ok(fluid_pressure);
    }

    /// get Re from pressure drop 
    #[inline] 
    pub fn get_reynolds_from_pressure_loss(&self,
        pressure_loss_input: Pressure,
        hydraulic_diameter: Length,
        fluid_density: MassDensity,
        fluid_viscosity:DynamicViscosity
    ) -> Result<Ratio,ThermalHydraulicsLibError>{

        if fluid_viscosity.value <= 0.0 {
            panic!("fluid Viscosity <= 0.0, nonphysical");
        }

        if hydraulic_diameter.value <= 0.0 {
            panic!("hydraulic Diameter <= 0.0, nonphysical");
        }

        if fluid_density.value <= 0.0 {
            panic!("fluidDensity <= 0.0, nonphysical");
        }

        // convert fluid pressure to bejan number 

        let bejan_input: Ratio = pressure_loss_input 
        * hydraulic_diameter
        * hydraulic_diameter
        * fluid_density
        / fluid_viscosity 
        / fluid_viscosity;

        // get the reynolds number 
        let reynolds_number_result = self.get_reynolds_number_from_bejan(
            bejan_input);

        reynolds_number_result
    }



}
