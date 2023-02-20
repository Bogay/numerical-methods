use bogay::{gram_schmidt, Matrix2D, Vec2};

fn main() {
    let a = Matrix2D::from_vec(
        Vec2::new(4, 5),
        vec![
            4., 2., 3., 0., -2., 3., -1., 1., 1., 3., -4., 2., 1., 0., 1., -1., 3., 1., 3., -2.,
        ],
    )
    .unwrap();
    println!("A =\n{a:>4}");
    let (q, r) = gram_schmidt(a, false);
    println!("Q =\n{q:>10.7}");
    println!("R =\n{r:>10.7}");

    let qr = q.mul(r).unwrap();
    println!("Q * R =\n{qr:>10.7}");
}
