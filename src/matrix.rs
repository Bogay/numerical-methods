use crate::vec2::{Square, Vec2};
use num_traits::One;
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub},
    str::FromStr,
};

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Matrix2D<T> {
    store: Vec<T>,
    size: Vec2,
}

impl<T> Matrix2D<T>
where
    T: Default + Clone,
{
    pub fn new(size: Vec2) -> Self {
        let store = vec![T::default(); size.x as usize * size.y as usize];
        Self { store, size }
    }
}

impl<T> Matrix2D<T>
where
    T: Clone,
{
    /// Create a matrix filled with given value
    pub fn fill(size: Vec2, value: T) -> Self {
        Self {
            size,
            store: vec![value; size.x as usize * size.y as usize],
        }
    }

    /// Try fill given area with given value, return error if the area is out of range
    pub fn try_fill(&mut self, anchor: Vec2, size: Vec2, value: T) -> Result<(), String> {
        let square = Square::new(anchor, size);
        // Check the fillin area is not out of range
        if square.row_iter().any(|pos| self.get(pos).is_none()) {
            return Err("Fill area out of range".to_string());
        }
        // Fillin
        for pos in square.row_iter() {
            *self.get_mut(pos).unwrap() = value.clone();
        }

        Ok(())
    }
}

impl<T> Matrix2D<T>
where
    T: Clone + Default + PartialEq,
{
    /// Try fill a given area without overwrite cells have already had default value
    pub fn try_fill_without_cover(
        &mut self,
        anchor: Vec2,
        size: Vec2,
        value: T,
    ) -> Result<(), String> {
        let square = Square::new(anchor, size);
        // Check there is not overwriting
        for pos in square.row_iter() {
            match self.get(pos) {
                Some(value) if value != &T::default() => {
                    return Err("Fill area covers non-default value".to_string())
                }
                None => return Err("Fill area out of range".to_string()),
                _ => {}
            }
        }
        // Fillin
        for pos in square.row_iter() {
            *self.get_mut(pos).unwrap() = value.clone();
        }

        Ok(())
    }
}

impl<T> Matrix2D<T>
where
    T: One + Default + Clone,
{
    pub fn identity(row: usize, col: usize) -> Self {
        let mut ret = Self::new(Vec2::new(col as i8, row as i8));
        for i in 0..::std::cmp::min(row, col) {
            *ret.get_mut(Vec2::new(i as i8, i as i8)).unwrap() = T::one();
        }
        ret
    }
}

impl<T> Matrix2D<T>
where
    T: Default + AddAssign + Clone + Mul<Output = T>,
{
    pub fn mul(&self, other: Self) -> Result<Self, String> {
        if self.size().x != other.size().y {
            return Err("Size of matrix does not match.".into());
        }

        let size = Vec2::new(other.size().x, self.size().y);
        let mut store = Vec::with_capacity(size.x as usize * size.y as usize);

        for r in 0..size.y {
            for c in 0..size.x {
                let mut acc = T::default();
                for i in 0..self.col() {
                    let i = i as i8;
                    let left = self[(i, r)].clone();
                    let right = other[(c, i)].clone();
                    let prod = (left * right).clone();
                    acc += prod;
                }
                store.push(acc);
            }
        }

        Self::from_vec(size, store)
    }
}

impl<T> Add for Matrix2D<T>
where
    T: Add<Output = T> + Clone,
{
    type Output = Matrix2D<T>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.size() != rhs.size() {
            panic!(
                "Different matrix size, except {}, got {}",
                self.size(),
                rhs.size()
            );
        }

        let mut store = Vec::with_capacity(self.col() * self.row());

        for i in 0..self.row() {
            let i = i as i8;
            for j in 0..self.col() {
                let j = j as i8;
                let a = self[(j, i)].clone();
                let b = rhs[(j, i)].clone();
                let r = a + b;
                store.push(r);
            }
        }

        Matrix2D::from_vec(self.size(), store).unwrap()
    }
}

impl<T> Sub for Matrix2D<T>
where
    T: Sub<Output = T> + Clone,
{
    type Output = Matrix2D<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        if self.size() != rhs.size() {
            panic!("Different matrix size");
        }

        let mut store = Vec::with_capacity(self.col() * self.row());

        for i in 0..self.row() {
            let i = i as i8;
            for j in 0..self.col() {
                let j = j as i8;
                let a = self[(j, i)].clone();
                let b = rhs[(j, i)].clone();
                let r = a - b;
                store.push(r);
            }
        }

        Matrix2D::from_vec(self.size(), store).unwrap()
    }
}

impl<T, S> MulAssign<S> for Matrix2D<T>
where
    T: MulAssign<S>,
    S: Copy,
{
    fn mul_assign(&mut self, rhs: S) {
        for ent in self.iter_mut() {
            *ent *= rhs;
        }
    }
}

impl<T> Matrix2D<T> {
    pub fn size(&self) -> Vec2 {
        self.size
    }

    pub fn row(&self) -> usize {
        self.size().y as usize
    }

    pub fn col(&self) -> usize {
        self.size().x as usize
    }

    /// Calculate the raw index of a Vec2 inside internal vec
    fn index(&self, pos: Vec2) -> usize {
        pos.y as usize * self.size.x as usize + pos.x as usize
    }

    /// Check whether the position is inside matrix
    fn is_inside(&self, pos: &Vec2) -> bool {
        pos.x >= 0 && pos.x < self.size.x && pos.y >= 0 && pos.y < self.size.y
    }

    /// Get a immutable reference of cell
    pub fn get(&self, pos: Vec2) -> Option<&T> {
        if !self.is_inside(&pos) {
            return None;
        }
        let idx = self.index(pos);
        self.store.get(idx)
    }

    /// Get a mutable reference of cell
    pub fn get_mut(&mut self, pos: Vec2) -> Option<&mut T> {
        if !self.is_inside(&pos) {
            return None;
        }
        let idx = self.index(pos);
        self.store.get_mut(idx)
    }

    /// Create matrix from given vector
    pub fn from_vec(size: Vec2, vec: Vec<T>) -> Result<Self, String> {
        let expect_size = size.x as usize * size.y as usize;
        if expect_size != vec.len() {
            return Err(format!(
                "Invalid vector size. expect {}, got {}",
                expect_size,
                vec.len()
            ));
        }

        Ok(Self { size, store: vec })
    }

    fn parse_size(line: &str) -> Result<Vec2, String> {
        let size = line.split_whitespace().collect::<Vec<_>>();
        if size.len() != 2 {
            return Err("First line should be the board row & column size".to_string());
        }
        let size = size
            .into_iter()
            .map(|s| {
                s.parse::<i8>()
                    .map_err(|e| format!("Failed to parse size: {}", e))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Vec2::new(size[1], size[0]))
    }

    /// Swap 2 rows R_i, R_j
    pub fn swap(&mut self, i: usize, j: usize) -> Result<(), String> {
        if i >= self.row() || j >= self.row() {
            return Err("Row index out of range".into());
        }

        for k in 0..self.col() {
            let a = self.index(Vec2::new(i as i8, k as i8));
            let b = self.index(Vec2::new(j as i8, k as i8));
            self.store.swap(a, b)
        }

        Ok(())
    }

    // Multiply any row by a constant
    pub fn mul_scalar<C>(&mut self, row: usize, c: C) -> Result<(), String>
    where
        T: MulAssign<C>,
        C: Clone,
    {
        if row >= self.row() {
            return Err("Row index out of range".into());
        }

        for i in 0..self.col() {
            *self.get_mut(Vec2::new(i as i8, row as i8)).unwrap() *= c.clone();
        }

        Ok(())
    }

    pub fn add_row<C>(&mut self, dst: usize, src: usize, c: C) -> Result<(), String>
    where
        T: AddAssign<T> + Mul<C, Output = T> + Clone,
        C: Clone,
    {
        if dst > self.row() || src > self.row() {
            return Err("Row index out of range".into());
        }

        for i in 0..self.col() {
            let i = i as i8;
            let v = self[(src as i8, i)].clone() * c.clone();
            self[(dst as i8, i)] += v;
        }

        Ok(())
    }

    pub fn transpose(self) -> Self {
        let mut rows = (0..self.col())
            .map(|_| Vec::with_capacity(self.row()))
            .collect::<Vec<_>>();
        let (row, col) = (self.row(), self.col());
        for (i, v) in self.store.into_iter().enumerate() {
            let i = i % rows.len();
            rows[i].push(v);
        }

        let store = rows.into_iter().flatten().collect();
        Self {
            store,
            size: Vec2::new(row as i8, col as i8),
        }
    }

    pub fn iter_row(&self) -> IterRow<T> {
        IterRow { mat: self, idx: 0 }
    }

    pub fn iter_col(&self) -> IterCol<T> {
        IterCol { mat: self, idx: 0 }
    }

    fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.store.iter_mut()
    }
}

pub struct IterRow<'a, T> {
    mat: &'a Matrix2D<T>,
    idx: usize,
}

impl<'a, T> Iterator for IterRow<'a, T> {
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.mat.row() {
            return None;
        }

        let mut v = Vec::with_capacity(self.mat.col());
        for i in 0..self.mat.col() {
            v.push(&self.mat[(i as i8, self.idx as i8)]);
        }

        self.idx += 1;
        Some(v)
    }
}

pub struct IterCol<'a, T> {
    mat: &'a Matrix2D<T>,
    idx: usize,
}

impl<'a, T> Iterator for IterCol<'a, T> {
    type Item = Vec<&'a T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.mat.col() {
            return None;
        }

        let mut v = Vec::with_capacity(self.mat.row());
        for i in 0..self.mat.row() {
            v.push(&self.mat[(self.idx as i8, i as i8)]);
        }

        self.idx += 1;
        Some(v)
    }
}

impl<T> FromStr for Matrix2D<T>
where
    T: FromStr,
    <T as FromStr>::Err: Debug,
{
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let mut input = input.lines();
        let line = input
            .next()
            .ok_or_else(|| "Missing first line".to_string())?;
        let size = Self::parse_size(line)?;

        if size.x <= 0 || size.y <= 0 {
            return Err("Either row or column size should >= 0".to_string());
        }

        let mut id_grid = Vec::with_capacity(size.x as usize * size.y as usize);
        for (row_i, line) in input.into_iter().take(size.y as usize).enumerate() {
            let row = line
                .split_whitespace()
                .map(|v| {
                    v.parse::<T>()
                        .map_err(|e| format!("Failed to parse block id: {:?}", e))
                })
                .collect::<Result<Vec<_>, _>>()?;
            if row.len() != size.x as usize {
                return Err(format!(
                    "Invalid line {}: expect {} block, got {}",
                    row_i,
                    size.x,
                    row.len(),
                ));
            }
            id_grid.extend(row);
        }
        Matrix2D::from_vec(size, id_grid)
    }
}

impl<T> Display for Matrix2D<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.row() {
            for j in 0..self.col() {
                let v = self.get(Vec2::new(j as i8, i as i8)).unwrap();
                f.pad(&format!("{} ", v))?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

impl<I0, I1, T> Index<(I0, I1)> for Matrix2D<T>
where
    I0: Into<i8>,
    I1: Into<i8>,
    // T: Clone,
{
    type Output = T;

    fn index(&self, (x, y): (I0, I1)) -> &Self::Output {
        self.get(Vec2::new(x.into(), y.into())).unwrap()
    }
}

impl<I0, I1, T> IndexMut<(I0, I1)> for Matrix2D<T>
where
    I0: Into<i8>,
    I1: Into<i8>,
    T: Clone,
{
    fn index_mut(&mut self, (x, y): (I0, I1)) -> &mut Self::Output {
        self.get_mut(Vec2::new(x.into(), y.into())).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eq() {
        let a = Matrix2D::fill(Vec2::new(3, 3), 5);
        let b = Matrix2D::fill(Vec2::new(3, 3), 5);

        assert_eq!(a, b);
    }

    #[test]
    fn test_ne() {
        let a = Matrix2D::fill(Vec2::new(3, 3), 2);
        let b = Matrix2D::fill(Vec2::new(3, 3), 5);

        assert_ne!(a, b);
    }

    #[test]
    fn test_get() -> Result<(), String> {
        let size = Vec2::new(3, 3);
        let mut v = vec![];
        for i in 0..9 {
            v.push(i);
        }
        let mat = Matrix2D::from_vec(size, v)?;

        for y in 0..3 {
            for x in 0..3 {
                assert_eq!(Some(&(y * 3 + x)), mat.get(Vec2::new(x, y)));
            }
        }

        Ok(())
    }

    #[test]
    fn test_get_out_of_range() {
        let mat = Matrix2D::fill(Vec2::new(2, 2), 7);

        assert_eq!(mat.get(Vec2::new(3, 1)), None);
        assert_eq!(mat.get(Vec2::new(1, 3)), None);
        assert_eq!(mat.get(Vec2::new(3, 3)), None);
    }

    #[test]
    fn test_add() {
        let a = Matrix2D::from_vec(Vec2::new(2, 2), vec![1, 2, 3, 4]).unwrap();
        let b = Matrix2D::from_vec(Vec2::new(2, 2), vec![5, 6, 7, 8]).unwrap();
        let c = a + b;

        let expected = Matrix2D::from_vec(Vec2::new(2, 2), vec![6, 8, 10, 12]).unwrap();
        assert_eq!(c, expected);
    }

    #[test]
    fn test_sub() {
        let a = Matrix2D::from_vec(Vec2::new(2, 2), vec![4, 3, 2, 1]).unwrap();
        let b = Matrix2D::from_vec(Vec2::new(2, 2), vec![5, 6, 7, 8]).unwrap();
        let c = b - a;

        let expected = Matrix2D::from_vec(Vec2::new(2, 2), vec![1, 3, 5, 7]).unwrap();
        assert_eq!(c, expected);
    }

    #[test]
    fn test_mul() {
        let a = Matrix2D::from_vec(Vec2::new(2, 2), vec![4, 3, 2, 1]).unwrap();
        let b = Matrix2D::from_vec(Vec2::new(2, 2), vec![5, 6, 7, 8]).unwrap();
        let c = a.mul(b).unwrap();

        let expected = Matrix2D::from_vec(Vec2::new(2, 2), vec![41, 48, 17, 20]).unwrap();
        assert_eq!(c, expected);

        let a = Matrix2D::from_vec(Vec2::new(1, 3), vec![1, 2, 2])
            .unwrap()
            .transpose();
        let b = Matrix2D::from_vec(Vec2::new(1, 3), vec![4, 3, 2]).unwrap();
        let c = a.mul(b).unwrap();

        let expected = Matrix2D::from_vec(Vec2::new(1, 1), vec![14]).unwrap();
        assert_eq!(c, expected);
    }

    #[test]
    fn test_tranpose() {
        let m = Matrix2D::from_vec(Vec2::new(2, 3), vec![1, 2, 3, 4, 5, 6])
            .unwrap()
            .transpose();
        let expected = Matrix2D::from_vec(Vec2::new(3, 2), vec![1, 3, 5, 2, 4, 6]).unwrap();
        assert_eq!(m, expected);
    }

    #[test]
    fn test_mul_scalar() {
        let mut m = Matrix2D::from_vec(Vec2::new(3, 2), vec![1, 2, 3, 4, 5, 6]).unwrap();
        m.mul_scalar(1, 3).unwrap();
        let expected = Matrix2D::from_vec(Vec2::new(3, 2), vec![1, 2, 3, 12, 15, 18]).unwrap();
        assert_eq!(m, expected);
    }

    #[test]
    fn test_iter_row() {
        let m = Matrix2D::from_vec(Vec2::new(3, 2), vec![1, 2, 3, 4, 5, 6]).unwrap();
        let mut it = m.iter_row();
        assert_eq!(it.next(), Some(vec![&1, &2, &3]));
        assert_eq!(it.next(), Some(vec![&4, &5, &6]));
        assert_eq!(it.next(), None);
    }

    #[test]
    fn test_iter_col() {
        let m = Matrix2D::from_vec(Vec2::new(3, 2), vec![1, 2, 3, 4, 5, 6]).unwrap();
        let mut it = m.iter_col();
        assert_eq!(it.next(), Some(vec![&1, &4]));
        assert_eq!(it.next(), Some(vec![&2, &5]));
        assert_eq!(it.next(), Some(vec![&3, &6]));
        assert_eq!(it.next(), None);
    }
}
