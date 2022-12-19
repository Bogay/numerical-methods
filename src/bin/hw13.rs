use num_traits::Pow;

fn main() {
    let f = |x: f64| -> f64 { x.pow(4) - 2. * x.pow(3) + 2. };
    let df = |x: f64| -> f64 { 4. * x.pow(3) - 6. * x.pow(2) };

    let eta = 0.0001;
    let mut x = 0.123;
    const EPSILON: f64 = 0.00000001;

    while df(x).abs() > EPSILON {
        x = x - df(x) * eta;
    }

    println!("x = {x}");
    println!("f(x) = {}", f(x));
}
