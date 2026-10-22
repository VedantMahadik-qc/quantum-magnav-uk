//! Synthetic anomaly fields for testing the filter before real data is involved.
//!
//! Every number used to build these fields is a **synthetic test value, not a physical
//! claim**. Amplitudes of a few hundred nT were chosen to be of the same order as the BGS
//! Great Britain survey (a sample spans about -990 to +1520 nT; see `docs/research-notes.md`).

use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

use crate::map::{Bounds, FieldMap, MapError, MapGrid};

/// One isotropic Gaussian anomaly: `A * exp(-r^2 / (2 sigma^2))`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianAnomaly {
    /// Centre, east (m).
    pub x: f64,
    /// Centre, north (m).
    pub y: f64,
    /// Peak value (nT); may be negative.
    pub amplitude_nt: f64,
    /// Width (m).
    pub sigma_m: f64,
}

impl GaussianAnomaly {
    /// Field of this anomaly alone at `(x, y)`.
    pub fn value_at(&self, x: f64, y: f64) -> f64 {
        let r2 = (x - self.x).powi(2) + (y - self.y).powi(2);
        self.amplitude_nt * (-r2 / (2.0 * self.sigma_m * self.sigma_m)).exp()
    }
}

/// Parameters for [`GaussianSumField::random`].
#[derive(Debug, Clone, PartialEq)]
pub struct RandomFieldParams {
    /// Region the field is defined over; anomaly centres are drawn uniformly inside it.
    pub bounds: Bounds,
    /// Number of anomalies.
    pub n_anomalies: usize,
    /// Range of |amplitude| (nT); each anomaly gets a random sign.
    pub amplitude_range_nt: (f64, f64),
    /// Range of widths (m).
    pub sigma_range_m: (f64, f64),
}

/// A sum of Gaussian anomalies, defined (non-`None`) only inside `bounds`.
///
/// Evaluating it directly gives an exact analytic field with no interpolation error;
/// [`GaussianSumField::to_grid`] samples it onto a [`MapGrid`].
#[derive(Debug, Clone, PartialEq)]
pub struct GaussianSumField {
    anomalies: Vec<GaussianAnomaly>,
    bounds: Bounds,
}

impl GaussianSumField {
    /// Builds a field from explicit anomalies.
    pub fn new(anomalies: Vec<GaussianAnomaly>, bounds: Bounds) -> Self {
        Self { anomalies, bounds }
    }

    /// Draws a random field; the same `seed` always gives the same field.
    pub fn random(params: &RandomFieldParams, seed: u64) -> Self {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        let b = params.bounds;
        let (a_lo, a_hi) = params.amplitude_range_nt;
        let (s_lo, s_hi) = params.sigma_range_m;
        let anomalies = (0..params.n_anomalies)
            .map(|_| {
                let sign = if rng.random::<bool>() { 1.0 } else { -1.0 };
                GaussianAnomaly {
                    x: rng.random_range(b.x_min..=b.x_max),
                    y: rng.random_range(b.y_min..=b.y_max),
                    amplitude_nt: sign * rng.random_range(a_lo..=a_hi),
                    sigma_m: rng.random_range(s_lo..=s_hi),
                }
            })
            .collect();
        Self::new(anomalies, b)
    }

    /// The individual anomalies.
    pub fn anomalies(&self) -> &[GaussianAnomaly] {
        &self.anomalies
    }

    /// Region where the field is defined.
    pub fn bounds(&self) -> Bounds {
        self.bounds
    }

    /// Exact field value, ignoring bounds.
    pub fn value_at(&self, x: f64, y: f64) -> f64 {
        self.anomalies.iter().map(|a| a.value_at(x, y)).sum()
    }

    /// Samples the field onto a grid with spacing `spacing_m` covering `bounds`.
    pub fn to_grid(&self, spacing_m: f64) -> Result<MapGrid, MapError> {
        let b = self.bounds;
        let nx = ((b.x_max - b.x_min) / spacing_m).floor() as usize + 1;
        let ny = ((b.y_max - b.y_min) / spacing_m).floor() as usize + 1;
        MapGrid::from_fn(b.x_min, b.y_min, spacing_m, spacing_m, nx, ny, |x, y| {
            self.value_at(x, y)
        })
    }
}

impl FieldMap for GaussianSumField {
    fn field_at(&self, x: f64, y: f64) -> Option<f64> {
        self.bounds.contains(x, y).then(|| self.value_at(x, y))
    }
}

/// A constant field: carries **no** position information. Used to prove the filter cannot
/// converge without map information.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlatField {
    /// The constant value (nT).
    pub value_nt: f64,
    /// Region where the field is defined.
    pub bounds: Bounds,
}

impl FieldMap for FlatField {
    fn field_at(&self, x: f64, y: f64) -> Option<f64> {
        self.bounds.contains(x, y).then_some(self.value_nt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const BOUNDS: Bounds = Bounds {
        x_min: 0.0,
        x_max: 10_000.0,
        y_min: 0.0,
        y_max: 8_000.0,
    };

    fn params() -> RandomFieldParams {
        RandomFieldParams {
            bounds: BOUNDS,
            n_anomalies: 10,
            amplitude_range_nt: (100.0, 500.0),
            sigma_range_m: (500.0, 1500.0),
        }
    }

    #[test]
    fn gaussian_has_analytic_peak_and_width() {
        let g = GaussianAnomaly {
            x: 100.0,
            y: 200.0,
            amplitude_nt: 300.0,
            sigma_m: 50.0,
        };
        assert_eq!(g.value_at(100.0, 200.0), 300.0);
        // At r = sigma the value is A * exp(-1/2), in any direction.
        assert_relative_eq!(g.value_at(150.0, 200.0), 300.0 * (-0.5f64).exp());
        assert_relative_eq!(g.value_at(100.0, 150.0), 300.0 * (-0.5f64).exp());
    }

    #[test]
    fn sum_field_is_sum_of_parts_and_none_outside() {
        let field = GaussianSumField::random(&params(), 7);
        let (x, y) = (1234.0, 5678.0);
        let manual: f64 = field.anomalies().iter().map(|a| a.value_at(x, y)).sum();
        assert_eq!(field.field_at(x, y), Some(manual));
        assert_eq!(field.field_at(-1.0, 100.0), None);
    }

    #[test]
    fn random_field_is_reproducible_and_within_ranges() {
        let a = GaussianSumField::random(&params(), 7);
        assert_eq!(a, GaussianSumField::random(&params(), 7));
        assert_ne!(a, GaussianSumField::random(&params(), 8));
        assert!(a.anomalies().iter().any(|g| g.amplitude_nt < 0.0));
        assert!(a.anomalies().iter().any(|g| g.amplitude_nt > 0.0));
        for g in a.anomalies() {
            assert!((100.0..=500.0).contains(&g.amplitude_nt.abs()));
            assert!((500.0..=1500.0).contains(&g.sigma_m));
            assert!(BOUNDS.contains(g.x, g.y));
        }
    }

    #[test]
    fn grid_matches_field_at_nodes() {
        let field = GaussianSumField::random(&params(), 7);
        let grid = field.to_grid(500.0).unwrap();
        assert_eq!((grid.nx(), grid.ny()), (21, 17));
        assert_eq!(grid.x_max(), BOUNDS.x_max);
        assert_eq!(grid.y_max(), BOUNDS.y_max);
        for (x, y) in [(0.0, 0.0), (500.0, 1500.0), (10_000.0, 8_000.0)] {
            assert_relative_eq!(
                grid.interpolate(x, y).unwrap(),
                field.value_at(x, y),
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn flat_field_is_constant_inside_bounds() {
        let flat = FlatField {
            value_nt: 50.0,
            bounds: BOUNDS,
        };
        assert_eq!(flat.field_at(1.0, 1.0), Some(50.0));
        assert_eq!(flat.field_at(9_999.0, 7_999.0), Some(50.0));
        assert_eq!(flat.field_at(10_001.0, 1.0), None);
    }
}
