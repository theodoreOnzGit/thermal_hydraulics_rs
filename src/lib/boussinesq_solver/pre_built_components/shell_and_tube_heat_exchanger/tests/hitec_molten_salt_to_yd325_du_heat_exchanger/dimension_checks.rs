


#[test]
pub fn check_shell_side_fluid_hydraulic_diameter(){
    use uom::si::f64::*;
    use uom::si::length::meter;
    use uom::si::ratio::ratio;

    // from Du's heat exchanger type, except we use one inner tube
    let tube_side_od = Length::new::<meter>(0.014);
    let number_of_tubes = 19_u32;
    let shell_side_id = Length::new::<meter>(0.1);

    let shell_side_fluid_hydraulic_diameter: Length = 
        (shell_side_id * shell_side_id - number_of_tubes as f64 *
         tube_side_od * tube_side_od)/
        (shell_side_id + number_of_tubes as f64 * tube_side_od);

    let hydraulic_diameter_by_l: Ratio = Ratio::new::<ratio>(0.009);
    let pipe_length = Length::new::<meter>(1.95);

    let given_de_by_du: Length = hydraulic_diameter_by_l * pipe_length;

    // max error is 3%
    approx::assert_relative_eq!(
        given_de_by_du.get::<meter>(),
        shell_side_fluid_hydraulic_diameter.get::<meter>(),
        max_relative = 0.03,
        );
}
