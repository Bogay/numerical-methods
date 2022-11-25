use bogay::{jacobi, Matrix2D, Vec2};

fn input_a() -> Matrix2D<f64> {
    let mut ret = Matrix2D::new(Vec2::new(10, 10));

    for i in 0..10 {
        for j in 0..10 {
            *ret.get_mut(Vec2::new(j, i)).unwrap() = if i == j {
                3.
            } else if i8::abs(i - j) == 1 {
                -1.
            } else {
                0.
            };
        }
    }

    ret
}

fn main() {
    let a = input_a();
    let b = Matrix2D::from_vec(
        Vec2::new(1, 10),
        vec![2., 1., 1., 1., 1., 1., 1., 1., 1., 2.],
    )
    .unwrap();
    let x = jacobi(a, b).unwrap();
    println!("{}", x);
}
