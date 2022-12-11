extern crate bogay;

use bogay::{plu_solve, Matrix2D, Vec2};
use num_traits::Pow;

// sqrt( (x - x1) ^ 2 + (y - y1) ^ 2 ) - (r - k)
trait Hw12Func: Fn(f64, f64, f64, f64, f64, f64) -> f64 {}

impl<T: Fn(f64, f64, f64, f64, f64, f64) -> f64> Hw12Func for T {}

struct Funcs {}

impl Funcs {
    fn f() -> impl Hw12Func {
        |x1, y1, r, x, y, k| {
            let dx: f64 = x - x1;
            let dy: f64 = y - y1;
            let z: f64 = dx.pow(2) + dy.pow(2);

            z.pow(0.5) - (r + k)
        }
    }

    fn dx() -> impl Hw12Func {
        |x1, y1, r, x, y, k| {
            let p: f64 = x - x1;
            let q: f64 = y - y1;

            let p = p * p;
            let q = q * q;

            assert_ne!(p + q, 0.);

            (x - x1) / (p + q).pow(0.5)
        }
    }

    fn dy() -> impl Hw12Func {
        |x1, y1, r, x, y, k| {
            let p: f64 = x - x1;
            let q: f64 = y - y1;

            let p = p * p;
            let q = q * q;

            assert_ne!(p + q, 0.);

            (y - y1) / (p + q).pow(0.5)
        }
    }

    fn dr() -> impl Hw12Func {
        |x1, y1, r, x, y, k| -1.
    }
}

struct MyFunc {
    x: f64,
    y: f64,
    r: f64,
    f: Box<dyn Hw12Func>,
}

impl MyFunc {
    pub fn new<T: Hw12Func + 'static>(x: f64, y: f64, r: f64, f: T) -> Self {
        Self {
            x,
            y,
            r,
            f: Box::new(f),
        }
    }

    pub fn eval(&self, x: f64, y: f64, k: f64) -> f64 {
        (self.f)(self.x, self.y, self.r, x, y, k)
    }
}

fn apply(f: &Matrix2D<MyFunc>, x: &Matrix2D<f64>) -> Matrix2D<f64> {
    let mut r = vec![];
    let (x, y, k) = (x[(0, 0)], x[(1, 0)], x[(2, 0)]);

    for row in f.iter_row() {
        for f in row {
            r.push(f.eval(x, y, k));
        }
    }

    Matrix2D::from_vec(f.size(), r).unwrap()
}

fn multivariate_newton(
    d_f: Matrix2D<MyFunc>,
    f: Matrix2D<MyFunc>,
    mut x: Matrix2D<f64>,
) -> Matrix2D<f64> {
    for _ in 0..1000 {
        let d_r = apply(&d_f, &x);
        let f_r = apply(&f, &x);

        let lx = d_r.clone().transpose().mul(d_r.clone()).unwrap();
        let mut rx = d_r.clone().transpose().mul(f_r.clone()).unwrap();

        for i in 0..rx.row() {
            for j in 0..rx.col() {
                *rx.get_mut(Vec2::new(j as i8, i as i8)).unwrap() *= -1.;
            }
        }

        let v = plu_solve(lx, rx).transpose();
        // let diff: f64 = v.iter_row().map(|r| r[0]).sum();
        // println!("diff = {diff}");
        x = x + v;
        // println!("x =\n{x}");
    }

    x
}

fn main() {
    let inputs = [
        [-1., 0., 1.], // (x, y, r)
        [1., 0.5, 0.5],
        [1., -0.5, 0.5],
        [0., 1., 0.5],
    ];
    let mut df_mat = vec![];
    let mut f_mat = vec![];
    for circle in inputs {
        let [x, y, r] = circle;
        df_mat.push(MyFunc::new(x, y, r, Funcs::dx()));
        df_mat.push(MyFunc::new(x, y, r, Funcs::dy()));
        df_mat.push(MyFunc::new(x, y, r, Funcs::dr()));
        f_mat.push(MyFunc::new(x, y, r, Funcs::f()));
    }
    let df_mat = Matrix2D::from_vec(Vec2::new(3, 4), df_mat).unwrap();
    let f_mat = Matrix2D::from_vec(Vec2::new(1, 4), f_mat).unwrap();
    let x = Matrix2D::from_vec(Vec2::new(3, 1), vec![0.5, 0.5, 0.1]).unwrap();

    let ans = multivariate_newton(df_mat, f_mat, x);

    let x = ans.get(Vec2::new(0, 0)).unwrap();
    let y = ans.get(Vec2::new(1, 0)).unwrap();
    let r = ans.get(Vec2::new(2, 0)).unwrap();

    println!("x = {x:<8.6}");
    println!("y = {y:<8.6}");
    println!("r = {r:<8.6}");
}
