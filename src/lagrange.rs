use num_traits::Float;

/// Lagrange interpolation
pub struct Lagrange<N: Float> {
    points: Vec<(N, N)>,
}

impl<N: Float> Lagrange<N> {
    pub fn new(points: Vec<(N, N)>) -> Self {
        Self { points }
    }

    pub fn eval(&self, in_x: N) -> N {
        let mut result = N::zero();

        for (i, (x, y)) in self.points.iter().enumerate() {
            let p = self
                .points
                .iter()
                .enumerate()
                .map(|(j, (_x, _))| if i != j { in_x - *_x } else { N::one() })
                .reduce(|a, b| a * b)
                .unwrap();
            let q = self
                .points
                .iter()
                .enumerate()
                .map(|(j, (_x, _))| if i != j { *x - *_x } else { N::one() })
                .reduce(|a, b| a * b)
                .unwrap();
            let d = *y * (p / q);
            result = result + d;
        }

        result
    }
}
