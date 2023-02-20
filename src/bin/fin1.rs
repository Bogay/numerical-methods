use bogay::{least_square_approximation, Matrix2D, Vec2};

fn main() {
    let a = Matrix2D::from_vec(Vec2::new(1, 3), vec![4., 7., 11.]).unwrap();
    let b = Matrix2D::from_vec(Vec2::new(1, 3), vec![3., 5., 8.]).unwrap();

    let x = least_square_approximation(a, b);
    println!("{x}");
}
