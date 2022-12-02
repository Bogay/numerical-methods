use bogay::{gram_schmidt, Matrix2D, Vec2};

fn main() {
    let a = Matrix2D::from_vec(
        Vec2::new(3, 5),
        vec![
            3., -1., 2., 4., 1., 0., -3., 2., 1., 1., 1., 5., -2., 0., 3.,
        ],
    )
    .unwrap();
    println!("A =\n{}", a);
    let (q, r) = gram_schmidt(a, false);
    println!("Q =\n{}", q);
    println!("R =\n{}", r);
}
