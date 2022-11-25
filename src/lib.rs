mod matrix;
mod vec2;

pub use matrix::Matrix2D;
pub use vec2::Vec2;

use num_traits::{Float, Inv};
use rand::{distributions::Standard, prelude::*};
use std::ops::{Add, AddAssign, Sub};

pub fn jacobi<'a, T>(a: Matrix2D<T>, b: Matrix2D<T>) -> Result<Matrix2D<T>, String>
where
    T: Clone + Inv<Output = T> + Default + AddAssign + Float,
    Standard: Distribution<T>,
    Matrix2D<T>: Add<Output = Matrix2D<T>> + Sub<Output = Matrix2D<T>>,
{
    if a.col() != b.row() {
        return Err("a.col() != b.row()".into());
    }

    let mut rng = thread_rng();
    let mut x = b.clone();
    for i in 0..x.row() {
        *x.get_mut(Vec2::new(0, i as i8)).unwrap() = T::from(rng.gen_range(0.0..16.0)).unwrap();
    }

    let mut d_inv = Matrix2D::new(a.size());
    for i in 0..a.col() {
        let i = Vec2::new(i as i8, i as i8);
        *d_inv.get_mut(i).unwrap() = a.get(i).unwrap().clone().inv();
    }

    let mut l = a.clone();
    let mut u = a.clone();
    for i in 0..a.row() {
        for j in 0..a.col() {
            if i <= j {
                *l.get_mut(Vec2::new(j as i8, i as i8)).unwrap() = T::zero();
            }
            if i >= j {
                *u.get_mut(Vec2::new(j as i8, i as i8)).unwrap() = T::zero();
            }
        }
    }
    let lu = l + u;

    let mut done = false;
    while !done {
        let _x = d_inv.mul(b.clone() - lu.mul(x.clone()).unwrap()).unwrap();
        done = (0..x.row()).all(|i| {
            let pos = Vec2::new(0, i as i8);
            let diff = x.get(pos).unwrap().clone() - _x.get(pos).unwrap().clone();
            diff == T::zero()
        });
        x = _x;
    }

    Ok(x)
}
