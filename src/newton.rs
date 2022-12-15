use num_traits::Float;

// Newton interpolating
pub struct Newton<N: Float> {
    points: Vec<(N, N)>,
}

impl<N: Float> Newton<N> {
    pub fn new(points: Vec<(N, N)>) -> Self {
        Self { points }
    }

    pub fn eval(&self, in_x: N) -> N {
        let mut fs = self.points.iter().map(|(_, y)| *y).collect::<Vec<_>>();
        let xs = self.points.iter().map(|(x, _)| *x).collect::<Vec<_>>();

        for i in 1..self.points.len() {
            for j in (i..self.points.len()).rev() {
                fs[j] = (fs[j] - fs[j - 1]) / (xs[j] - xs[j - i]);
            }
        }

        let mut ans = N::zero();
        let mut x = N::one();
        for (f, xp) in fs.iter().zip(&xs) {
            ans = ans + (*f * x);
            x = x * (in_x - *xp);
        }

        ans
    }
}
