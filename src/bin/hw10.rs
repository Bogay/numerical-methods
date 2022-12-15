extern crate bogay;
use bogay::{least_square_approximation, Matrix2D, Vec2};

// Least square approximation
fn main() {
    // f(x) = a + bx + cx^2 + dx^4

    // x = range(0, 3, 0.25)
    let xs = (0..=12).map(|x| 0.25 * x as f64).collect::<Vec<_>>();
    // Build mat A
    let a_rows = xs
        .iter()
        .flat_map(|&x| [1., x, x * x, x * x * x])
        .collect::<Vec<_>>();
    let a = Matrix2D::from_vec(Vec2::new(4, xs.len() as i8), a_rows).unwrap();

    let ys = [
        6.3806, 7.1338, 9.1662, 11.5545, 15.6414, 22.7371, 32.0696, 47.0756, 73.1596, 111.4684,
        175.9895, 278.5550, 446.4441,
    ];
    let b = Matrix2D::from_vec(Vec2::new(1, ys.len() as i8), ys.to_vec()).unwrap();

    let x = least_square_approximation(a, b);
    let symbols = ['a', 'b', 'c', 'd'];
    for (i, sym) in symbols.into_iter().enumerate() {
        println!("{} = {}", sym, x.get(Vec2::new(0, i as i8)).unwrap());
    }
}
