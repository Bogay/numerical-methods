struct La {
    points: Vec<(f64, f64)>,
}

impl La {
    fn new(data_points: &[(f64, f64)]) -> Self {
        Self {
            points: data_points.into(),
        }
    }

    fn eval(&self, in_x: f64) -> f64 {
        let mut result = 0.;

        for (i, (x, y)) in self.points.iter().enumerate() {
            let p = self
                .points
                .iter()
                .enumerate()
                .map(|(j, (_x, _y))| if i != j { in_x - _x } else { 1. })
                .reduce(|a, b| a * b)
                .unwrap();
            let q = self
                .points
                .iter()
                .enumerate()
                .map(|(j, (_x, _y))| if i != j { x - _x } else { 1. })
                .reduce(|a, b| a * b)
                .unwrap();
            result += y * p / q;
        }

        result
    }
}

/// Lagrange interpolation
fn main() {
    let data_points = vec![
        (1.00, 0.1924),
        (1.05, 0.2414),
        (1.10, 0.2933),
        (1.15, 0.3492),
    ];
    let v = La::new(&data_points).eval(1.09);
    println!("{}", v);
}
