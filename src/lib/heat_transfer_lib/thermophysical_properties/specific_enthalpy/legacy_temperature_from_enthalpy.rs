/// returns temperature of stainless steel 304L 
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// The algorithm is to make a spline of enthalpy and temperature
/// and then use the enthalpy to obtain a temperature
///
/// Note: h in this function represents specific enthalpy
#[inline]
fn steel_304_l_spline_temp_attempt_2_from_specific_enthalpy(
    h_steel: AvailableEnergy) -> ThermodynamicTemperature {

    // the idea is basically to evaluate enthalpy at the 
    // following temperatures
    let temperature_values_kelvin: Vec<f64>
    = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);

    // and then use that to formulate a spline,
    // with the spline, i'll evaluate enthalpy from temperature
    // within pretty much one iteration. However, it is spline 
    // construction which may take a little long. 
    //
    // However, the number of iterations per calculation is fixed
    //
    // I won't optimise it now just yet

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let mut enthalpy_vector = vec![0.0; temperature_vec_len];

    for index_i in 0..temperature_vec_len {

        // first, evaluate the enthalpy at temperature values 
        let temperature_value = temperature_values_kelvin[index_i];

        //next let's evaluate the specific enthalpy of steel 
        let steel = Material::Solid(SteelSS304L);
        let steel_temp = ThermodynamicTemperature::new::<kelvin>(
            temperature_value);
        let pressure = Pressure::new::<atmosphere>(1.0);

        let steel_enthalpy_result = specific_enthalpy(steel, 
            steel_temp, pressure);

        let steel_enthalpy_value = match steel_enthalpy_result {
            Ok(steel_enthalpy) => steel_enthalpy.value,
            // i can of course unwrap the result,
            // but i want to leave it more explicit in case 
            // i wish to manually handle the error
            Err(error_msg) => panic!("{}",error_msg),
        };

        // once i evalute the enthalpy value, pass it on to the vector

        enthalpy_vector[index_i] = steel_enthalpy_value;

    }


    // now I have my enthalpy vector, i can do an inverted spline 
    // to have enthalpy given in as an input, and temperature received
    // as an output

    let enthalpy_to_temperature_spline = 
    CubicSpline::from_nodes(&enthalpy_vector,
    &temperature_values_kelvin);

    // now let's get our enthalpy in joules_per_kg
    let h_steel_joules_per_kg = h_steel.get::<joule_per_kilogram>();

    let temperature_from_enthalpy_kelvin = 
    enthalpy_to_temperature_spline.eval(h_steel_joules_per_kg);

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}

#[test]
pub fn steel_temperature_from_enthalpy_test_spline_2(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let steel = Material::Solid(SteelSS304L);
    let steel_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        steel,steel_temp,pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    steel_304_l_spline_temp_attempt_2_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically off by less than 0.5K, which is 
    // within measurement error!
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=0.5);


}

/// Old attempt of spline,
/// attempt 1
///
/// returns temperature of stainless steel 304L 
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// The algorithm is to make a spline of enthalpy and temperature
///
/// Note: h in this function represents specific enthalpy
///
/// I suspected that the integral function may have messed up the 
/// spline, thus we have unsatisfactory behaviour
#[inline]
fn steel_304_l_spline_temp_attempt_1_from_specific_enthalpy(
    h_steel: AvailableEnergy) -> ThermodynamicTemperature {

    let temperature_values_kelvin: Vec<f64>
    = c!(250.0, 300.0, 350.0, 
        400.0, 450.0, 500.0, 700.0, 1000.0);

    let specific_heat_capacity_values_joule_per_kilogram_kelvin = c!(443.3375,
        457.0361, 469.4894, 480.6974, 490.66, 500.6227, 526.7746,
        551.6812);

    let h_spline_integral = 
    CubicSpline::from_nodes(&temperature_values_kelvin, 
        &specific_heat_capacity_values_joule_per_kilogram_kelvin)
        .integral();

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let zero_celsius_vec = vec![273.15; temperature_vec_len];

    // now the integrals are indefinite, hence I need to ensure 
    // the enthalpies are taken with reference to the common 
    // reference point of 0 degrees C or 273.15 K
    let integral_spline_enthalpy_vector = 
    h_spline_integral.eval_vec(&temperature_values_kelvin);
    let reference_enthalpy_vector =
    h_spline_integral.eval_vec(&zero_celsius_vec);

    // this is a concise but rather complicated way of 
    // subtraction
    let enthalpy_vector: Vec<f64> = integral_spline_enthalpy_vector
        .iter()
        .zip(reference_enthalpy_vector)
        .map(|(elem_a, elem_b)| elem_a - elem_b)
        .collect();
    
    // now I have my enthalpy vector, i can do an inverted spline 
    // to have enthalpy given in as an input, and temperature received
    // as an output

    let enthalpy_to_temperature_spline = 
    CubicSpline::from_nodes(&enthalpy_vector,
    &temperature_values_kelvin);

    // now let's get our enthalpy in joules_per_kg
    let h_steel_joules_per_kg = h_steel.get::<joule_per_kilogram>();

    let temperature_from_enthalpy_kelvin = 
    enthalpy_to_temperature_spline.eval(h_steel_joules_per_kg);

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}

#[test]
pub fn steel_temperature_from_enthalpy_test_spline_1(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let steel = Material::Solid(SteelSS304L);
    let steel_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        steel,steel_temp,pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    steel_304_l_spline_temp_attempt_1_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically off by about 21K,
    // the temperature was i think 354K
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=21.0);

    approx::assert_abs_diff_ne!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=0.0);

    // looks like the spline is good as an initial guess, but not 
    // really for finding the answer

}

/// returns temperature of copper
/// cited from:
/// Zou, L., Hu, R., & Charpentier, A. (2019). SAM code 
/// validation using the compact integral effects test (CIET) experimental 
/// data (No. ANL/NSE-19/11). Argonne National 
/// Lab.(ANL), Argonne, IL (United States).
///
/// The algorithm is to make a spline of enthalpy and temperature
/// and then use the enthalpy to obtain a temperature
///
/// Note: h in this function represents specific enthalpy
#[inline]
fn copper_spline_temp_attempt_1_from_specific_enthalpy(
    h_copper: AvailableEnergy) -> ThermodynamicTemperature {

    // the idea is basically to evaluate enthalpy at the 
    // following temperatures
    let temperature_values_kelvin: Vec<f64>
    = c!(200.0 ,250.0, 300.0, 350.0, 
        400.0, 500.0, 1000.0);

    // and then use that to formulate a spline,
    // with the spline, i'll evaluate enthalpy from temperature
    // within pretty much one iteration. However, it is spline 
    // construction which may take a little long. 
    //
    // However, the number of iterations per calculation is fixed
    //
    // I won't optimise it now just yet

    let temperature_vec_len = 
    temperature_values_kelvin.len();

    let mut enthalpy_vector = vec![0.0; temperature_vec_len];

    for index_i in 0..temperature_vec_len {

        // first, evaluate the enthalpy at temperature values 
        let temperature_value = temperature_values_kelvin[index_i];

        //next let's evaluate the specific enthalpy of copper 
        let copper = Material::Solid(Copper);
        let copper_temp = ThermodynamicTemperature::new::<kelvin>(
            temperature_value);
        let pressure = Pressure::new::<atmosphere>(1.0);

        let copper_enthalpy_result = specific_enthalpy(copper, 
            copper_temp, pressure);

        let copper_enthalpy_value = match copper_enthalpy_result {
            Ok(copper_enthalpy) => copper_enthalpy.value,
            // i can of course unwrap the result,
            // but i want to leave it more explicit in case 
            // i wish to manually handle the error
            Err(error_msg) => panic!("{}",error_msg),
        };

        // once i evalute the enthalpy value, pass it on to the vector

        enthalpy_vector[index_i] = copper_enthalpy_value;

    }


    // now I have my enthalpy vector, i can do an inverted spline 
    // to have enthalpy given in as an input, and temperature received
    // as an output

    let enthalpy_to_temperature_spline = 
    CubicSpline::from_nodes(&enthalpy_vector,
    &temperature_values_kelvin);

    // now let's get our enthalpy in joules_per_kg
    let h_copper_joules_per_kg = h_copper.get::<joule_per_kilogram>();

    let temperature_from_enthalpy_kelvin = 
    enthalpy_to_temperature_spline.eval(h_copper_joules_per_kg);

    // return temperature
    ThermodynamicTemperature::new::<kelvin>(
        temperature_from_enthalpy_kelvin)

}

#[test]
pub fn copper_temperature_from_enthalpy_test_spline_1(){
    // we'll test temperature at 375K 
    // we should get an enthalpy from the spline 
    // for zweibaum's paper 

    let copper = Material::Solid(Copper);
    let copper_temp = ThermodynamicTemperature::new::<kelvin>(375.0);
    let pressure = Pressure::new::<atmosphere>(1.0);

    let enthalpy_spline_zweibaum_375k = specific_enthalpy(
        copper,copper_temp,pressure).unwrap();

    // now we have an enthalpy, let's check the temperature 

    let temperature_from_enthalpy_test = 
    copper_spline_temp_attempt_1_from_specific_enthalpy(
        enthalpy_spline_zweibaum_375k);

    // we are basically by about 5K, which is 
    // not within measurement error, probably have to do more work
    // what this means is that accuracy is sacrificed
    // for speed, sometimes too much accuracy
    //
    // for enthalpy, we probably want to have it as accurate 
    // as possible so that energy doesn't appear from nowhere 
    // and disappear from calculation
    //
    // I would note though, that the spline method does 
    // give a pretty good initial guess of where the temperature 
    // ought to be, so perhaps the iterative method can be used 
    // for the last few iterations to convergence
    // we could use brent dekker method
    approx::assert_abs_diff_eq!(
        temperature_from_enthalpy_test.get::<kelvin>(),
        375.0,
        epsilon=5.0);


}
