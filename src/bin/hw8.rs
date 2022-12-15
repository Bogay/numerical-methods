use bogay::Lagrange;

fn main() {
    let data_points = vec![
        (1., 10.),
        (2., 10.),
        (3., 10.),
        (4., 10.),
        (5., 10.),
        (6., 15.),
    ];
    let v = Lagrange::new(data_points).eval(7.);
    println!("{}", v);
}
