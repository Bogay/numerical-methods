use num_traits::{Float, Pow};
use rand::prelude::*;

#[derive(Debug, Clone)]
struct Eta<T: Float> {
    val: T,
    min: T,
    decay: T,
}

impl<T: Float> Eta<T> {
    pub fn new(val: T, min: T, decay: T) -> Self {
        assert!(val >= min);
        assert!(decay > T::zero());

        return Self { val, min, decay };
    }

    pub fn update(&mut self) {
        self.val = (self.val - self.decay).max(self.min);
    }

    pub fn val(&self) -> T {
        self.val
    }
}

fn f(x: f64) -> f64 {
    x.pow(4) + 3. * (x - 2.).pow(3) - 15. * x.pow(2) + 1.
}

fn df(x: f64) -> f64 {
    4. * x.pow(3) + 9. * x.pow(2) - 66. * x + 36.
}

fn main() {
    const EPSILON: f64 = 0.000000001;
    let approx = |mut eta: Eta<f64>, mut x: f64| {
        let mut t = 1000000;
        while df(x).abs() > EPSILON {
            x = x - df(x) * eta.val();
            t -= 1;
            if t == 0 {
                break;
            }
            eta.update();
        }
        x
    };

    let eta = [
        Eta::new(0.0001, 0.000001, 0.000001),
        Eta::new(0.00001, 0.00000001, 0.000000001),
        Eta::new(0.001, 0.0000001, 0.000001),
        Eta::new(0.0001, 0.00001, 0.0000001),
    ];

    let mut rng = thread_rng();
    let x = (0..15000)
        .map(|_| {
            approx(
                eta.choose(&mut rng).unwrap().clone(),
                rng.gen_range((-500.)..500.),
            )
        })
        // filter NaN and Inf
        .filter(|x| x.is_finite())
        // can't use min* due to f64 does not impl Ord
        .reduce(|xp, mut x| {
            if f(x) < f(xp) {
                x = xp;
            }
            x
        })
        .unwrap();

    println!("x = {x}");
    println!("f(x) = {}", f(x));
}
