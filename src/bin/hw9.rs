struct Newton {
    points: Vec<(f64, f64)>,
}

impl Newton {
    fn new(data_points: &[(f64, f64)]) -> Self {
        Self {
            points: data_points.into(),
        }
    }

    fn eval(&self, in_x: f64) -> f64 {
        let mut fs = self.points.iter().map(|(_, y)| *y).collect::<Vec<_>>();
        let xs = self.points.iter().map(|(x, _)| *x).collect::<Vec<_>>();

        for i in 1..self.points.len() {
            for j in (i..self.points.len()).rev() {
                fs[j] = (fs[j] - fs[j - 1]) / (xs[j] - xs[j - i]);
            }
        }

        let mut ans = 0.;
        let mut x = 1.;
        for (f, xp) in fs.iter().zip(&xs) {
            ans += f * x;
            x *= in_x - xp;
        }

        ans
    }
}

/// Use Newton’s divided differences to find the degree 4 interpolating polynomial for the data
fn main() {
    let data_points = vec![
        (0.6, 1.433329),
        (0.7, 1.632316),
        (0.8, 1.896481),
        (0.9, 2.247908),
        (1.0, 2.718282),
    ];
    for x in [0.82, 0.98] {
        let y = Newton::new(&data_points).eval(x);
        println!("P({}) = {}", x, y);
    }
}
