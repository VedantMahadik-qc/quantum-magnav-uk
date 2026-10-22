//! A tiny 2-D vector type for positions (m) and velocities (m/s).
//!
//! Implementing the `std::ops` traits (`Add`, `Sub`, `Mul`) is how Rust does operator
//! overloading: after these impls, `a + b * 2.0` works on `Vec2` values.

use std::ops::{Add, AddAssign, Mul, Sub};

/// A 2-D vector: `x` = east, `y` = north.
///
/// `Copy` means it is duplicated on assignment like an `f64`, rather than moved, which is
/// cheap here (16 bytes) and avoids ownership bookkeeping for small values.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec2 {
    /// East component.
    pub x: f64,
    /// North component.
    pub y: f64,
}

impl Vec2 {
    /// Creates a vector from its components.
    pub const fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// Euclidean length.
    pub fn norm(self) -> f64 {
        self.x.hypot(self.y)
    }
}

impl Add for Vec2 {
    type Output = Vec2;
    fn add(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl AddAssign for Vec2 {
    fn add_assign(&mut self, rhs: Vec2) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Sub for Vec2 {
    type Output = Vec2;
    fn sub(self, rhs: Vec2) -> Vec2 {
        Vec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Vec2 {
    type Output = Vec2;
    fn mul(self, k: f64) -> Vec2 {
        Vec2::new(self.x * k, self.y * k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_and_norm() {
        let a = Vec2::new(3.0, 4.0);
        let b = Vec2::new(1.0, -2.0);
        assert_eq!(a + b, Vec2::new(4.0, 2.0));
        assert_eq!(a - b, Vec2::new(2.0, 6.0));
        assert_eq!(a * 2.0, Vec2::new(6.0, 8.0));
        assert_eq!(a.norm(), 5.0);
        let mut c = a;
        c += b;
        assert_eq!(c, a + b);
    }
}
