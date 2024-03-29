mod lagrange;
mod matrix;
mod newton;
mod vec2;

pub use lagrange::Lagrange;
pub use matrix::Matrix2D;
pub use newton::Newton;
pub use vec2::Vec2;

use num_traits::{Float, Inv};
use rand::{distributions::Standard, prelude::*};
use std::{
    fmt::Display,
    ops::{Add, AddAssign, Sub},
};

pub fn jacobi<T>(a: Matrix2D<T>, b: Matrix2D<T>) -> Result<Matrix2D<T>, String>
where
    T: Inv<Output = T> + Default + AddAssign + Float + Display,
    Standard: Distribution<T>,
    Matrix2D<T>: Add<Output = Matrix2D<T>> + Sub<Output = Matrix2D<T>>,
{
    if a.col() != b.row() {
        return Err("a.col() != b.row()".into());
    }

    let mut rng = thread_rng();
    let mut x = b.clone();
    for i in 0..x.row() {
        *x.get_mut(Vec2::new(0, i as i8)).unwrap() = T::from(rng.gen_range(-5.0..=5.0)).unwrap();
    }

    let mut d_inv = Matrix2D::new(a.size());
    for i in 0..a.col() {
        let i = i as i8;
        d_inv[(i, i)] = a[(i, i)].inv();
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

    loop {
        let _x = d_inv.mul(b.clone() - lu.mul(x.clone()).unwrap()).unwrap();

        let done = (0..x.row()).all(|i| {
            let diff = x[(0, i as i8)] - _x[(0, i as i8)];
            diff == T::zero()
        });
        if done {
            break;
        }

        x = _x;
    }

    Ok(x)
}

pub fn plu_solve(a: Matrix2D<f64>, mut b: Matrix2D<f64>) -> Matrix2D<f64> {
    assert_eq!(a.row(), b.row());
    assert_eq!(b.col(), 1);

    let mut u = a.clone();
    let mut l = Matrix2D::new(a.size());
    let mut p = Matrix2D::<f64>::identity(a.row(), a.col());

    for c in 0..u.col() {
        let c = c as i8;
        if u[(c, c)] == 0. {
            let mut k = c + 1;
            let row = u.row() as i8;
            while k < row {
                if u.get(Vec2::new(c, k)).unwrap() != &0. {
                    u.swap(k as usize, c as usize).unwrap();
                    p.swap(k as usize, c as usize).unwrap();
                    l.swap(k as usize, c as usize).unwrap();
                    b.swap(k as usize, c as usize).unwrap();
                    break;
                }
                k += 1;
            }
            assert_ne!(k, row);
        }

        let d = *u.get(Vec2::new(c, c)).unwrap();
        assert_ne!(d, 0.);
        for r in c + 1..u.row() as i8 {
            let s = -u[(c, r)] / d;
            u.add_row(r as usize, c as usize, s).unwrap();
            l[(c, r)] = -s;
        }
    }

    for i in 0..l.row().min(l.col()) {
        l[(i as i8, i as i8)] = 1.;
    }

    let mut x = b.clone();
    // solve Lc = B
    for r in 0..l.row() {
        let r = r as i8;
        for c in 0..r {
            let c = c as i8;
            let v = x[(0, c)] * l[(c, r)];
            x[(0, r)] -= v;
        }
        x[(0, r)] /= l[(r, r)];
    }
    // solve Ux = c
    for r in (0..u.row()).rev() {
        let r = r as i8;
        for c in r + 1..u.row() as i8 {
            let v = x[(0, c)] * u[(c, r)];
            x[(0, r)] -= v;
        }

        x[(0, r)] /= u[(r, r)];
    }

    x
}

/// Apply Gram-Schmidt orthogonalization to find reduced QR factorization on `a`
pub fn gram_schmidt(a: Matrix2D<f64>, is_full_qr: bool) -> (Matrix2D<f64>, Matrix2D<f64>) {
    let mut a_cols: Vec<_> = a
        .iter_col()
        .map(|col| {
            Matrix2D::from_vec(
                Vec2::new(1, a.row() as i8),
                col.into_iter().copied().collect(),
            )
            .unwrap()
        })
        .collect();

    if is_full_qr {
        assert!(a.row() > a.col());
        let sz = a_cols[0].size();
        for _ in 0..(a.row() - a.col()) {
            a_cols.push(Matrix2D::fill(sz, 1.));
        }
    }

    let mut qs: Vec<Matrix2D<_>> = vec![];
    let mut r = Matrix2D::<f64>::new(Vec2::new(a.col() as i8, a_cols.len() as i8));
    for (j, aj) in a_cols.iter().enumerate() {
        let mut y = aj.clone();
        for (i, q) in qs.iter().enumerate() {
            let mut q = q.clone();
            let r_s = *q
                .clone()
                .transpose()
                .mul(aj.clone())
                .unwrap()
                .get(Vec2::new(0, 0))
                .unwrap();

            // It may out of range if is_full_qr = true
            if let Some(ent) = r.get_mut(Vec2::new(j as i8, i as i8)) {
                *ent = r_s;
            }

            q *= r_s;
            y = y - q;
        }

        let y_len = y.iter_row().map(|v| v[0] * v[0]).sum::<f64>().sqrt();
        assert_ne!(y_len, 0.);
        // It may out of range if is_full_qr = true
        if let Some(ent) = r.get_mut(Vec2::new(j as i8, j as i8)) {
            *ent = y_len;
        }
        y *= 1. / y_len;
        qs.push(y);
    }

    let mut q = Matrix2D::new(Vec2::new(qs.len() as i8, qs[0].row() as i8));
    for (i, iq) in qs.iter().enumerate() {
        for j in 0..iq.row() {
            let v = iq.get(Vec2::new(0, j as i8)).unwrap();
            *q.get_mut(Vec2::new(i as i8, j as i8)).unwrap() = *v;
        }
    }

    (q, r)
}

// Find least square solution x' for in consisten system Ax = b
pub fn least_square_approximation(a: Matrix2D<f64>, b: Matrix2D<f64>) -> Matrix2D<f64> {
    // Calculate A^T * A
    let at = a.clone().transpose();
    let ata = at.mul(a).unwrap();
    // Calculat A^T * b
    let bp = at.mul(b).unwrap();

    plu_solve(ata, bp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolation() {
        const EPSILON: f64 = 0.000001;
        let points = vec![
            (0.6, 1.433329),
            (0.7, 1.632316),
            (0.8, 1.896481),
            (0.9, 2.247908),
            (1.0, 2.718282),
        ];
        let ans = (0.82, 1.9589097744);

        let newton_eval = Newton::new(points.clone()).eval(ans.0);
        let lagrange_eval = Lagrange::new(points.clone()).eval(ans.0);

        assert!((ans.1 - newton_eval).abs() < EPSILON);
        assert!((ans.1 - lagrange_eval).abs() < EPSILON);
    }

    #[test]
    fn test_jacobi() {
        let a = Matrix2D::from_vec(Vec2::new(2, 2), vec![4., 2., 3., 4.]).unwrap();
        let b = Matrix2D::from_vec(Vec2::new(1, 2), vec![8., 11.]).unwrap();
        let x = jacobi(a.clone(), b.clone()).unwrap();
        let bp = a.mul(x.clone()).unwrap();

        assert_eq!(b, bp);
    }

    #[test]
    fn test_plu_solve() {
        let a = Matrix2D::from_vec(Vec2::new(2, 2), vec![3., 7., 2., 1.]).unwrap();
        let b = Matrix2D::from_vec(Vec2::new(1, 2), vec![17., 4.]).unwrap();
        let x = plu_solve(a.clone(), b.clone());
        let bp = a.mul(x.clone()).unwrap();

        assert_eq!(b, bp);
    }
}
