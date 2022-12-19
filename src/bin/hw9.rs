use bogay::Newton;

/// Use Newton’s divided differences to find the degree 4 interpolating polynomial for the data
fn main() {
    #[allow(clippy::approx_constant)]
    let data_points = vec![
        (0.6, 1.433329),
        (0.7, 1.632316),
        (0.8, 1.896481),
        (0.9, 2.247908),
        (1.0, 2.718282), // f64::consts::E
    ];
    for x in [0.82, 0.98] {
        let y = Newton::new(data_points.clone()).eval(x);
        println!("P({x}) = {y}");
    }
}
