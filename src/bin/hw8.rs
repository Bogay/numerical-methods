use bogay::Lagrange;

fn main() {
    let data_points = vec![
        (1.00, 0.1924),
        (1.05, 0.2414),
        (1.10, 0.2933),
        (1.15, 0.3492),
    ];
    let v = Lagrange::new(data_points).eval(1.09);
    println!("{}", v);
}
