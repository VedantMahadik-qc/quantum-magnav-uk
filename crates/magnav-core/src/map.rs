//! Gridded magnetic anomaly maps and bilinear interpolation.
//!
//! A [`MapGrid`] stores field values on a regular grid. Grid conventions:
//!
//! - `values` has shape `(ny, nx)`: **row index = y (north), column index = x (east)**.
//! - **Row 0 is the southernmost row.** `y` increases with row index. (Many image-style numpy
//!   arrays put north at row 0; such arrays must be flipped before use.)
//! - `(x0, y0)` is the position of the node `values[[0, 0]]`, i.e. the south-west corner node.
//!   Node `(row, col)` sits at `(x0 + col * dx, y0 + row * dy)`.
//! - `NaN` values mean "no data" (e.g. offshore gaps in a real survey).
//!
//! ## File format
//! A map on disk is two files with the same stem: `foo.npy` (a 2-D float64 numpy array in the
//! orientation above) and a sidecar `foo.json`:
//!
//! ```json
//! { "x0": 400000.0, "y0": 300000.0, "dx": 1000.0, "dy": 1000.0,
//!   "crs": "EPSG:27700", "units": "nT", "altitude_m": 305.0, "source": "BGS aeromag" }
//! ```
//!
//! `altitude_m` and `source` are optional. `units` must be `"nT"`.

use std::path::{Path, PathBuf};

use ndarray::Array2;
use ndarray_npy::{ReadNpyError, WriteNpyError, read_npy, write_npy};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Anything that can report a magnetic anomaly value at a position.
///
/// The particle filter is generic over this trait, so it runs unchanged on a gridded map,
/// an analytic test field or a flat (information-free) field.
pub trait FieldMap {
    /// Total-field anomaly in nT at `(x, y)` metres, or `None` where the map has no data.
    fn field_at(&self, x: f64, y: f64) -> Option<f64>;
}

/// Errors from building or loading a [`MapGrid`].
#[derive(Debug, Error)]
pub enum MapError {
    /// Grid spacing must be finite and strictly positive.
    #[error("grid spacing must be finite and positive, got dx = {dx}, dy = {dy}")]
    InvalidSpacing {
        /// Spacing in x (m).
        dx: f64,
        /// Spacing in y (m).
        dy: f64,
    },
    /// Bilinear interpolation needs at least a 2 x 2 grid.
    #[error("grid must be at least 2 x 2 for interpolation, got nx = {nx}, ny = {ny}")]
    TooSmall {
        /// Number of columns.
        nx: usize,
        /// Number of rows.
        ny: usize,
    },
    /// The grid origin must be a finite coordinate.
    #[error("grid origin must be finite, got x0 = {x0}, y0 = {y0}")]
    InvalidOrigin {
        /// Origin x (m).
        x0: f64,
        /// Origin y (m).
        y0: f64,
    },
    /// Reading or writing a file failed.
    #[error("I/O error on {path}: {source}")]
    Io {
        /// File that failed.
        path: PathBuf,
        /// Underlying OS error.
        source: std::io::Error,
    },
    /// The `.npy` file could not be read as a 2-D float64 array.
    #[error("cannot read {path} as a 2-D float64 .npy array: {source}")]
    ReadNpy {
        /// File that failed.
        path: PathBuf,
        /// Underlying parse error.
        source: ReadNpyError,
    },
    /// The `.npy` file could not be written.
    #[error("cannot write {path}: {source}")]
    WriteNpy {
        /// File that failed.
        path: PathBuf,
        /// Underlying error.
        source: WriteNpyError,
    },
    /// The JSON sidecar is malformed or missing fields.
    #[error("invalid map metadata in {path}: {source}")]
    Json {
        /// File that failed.
        path: PathBuf,
        /// Underlying parse error.
        source: serde_json::Error,
    },
    /// Only nanotesla maps are supported, to avoid silent unit mix-ups.
    #[error("map units must be \"nT\", got {0:?}")]
    WrongUnits(String),
}

/// Descriptive metadata stored alongside a map (everything except the grid geometry).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapMetadata {
    /// Coordinate reference system, e.g. `"EPSG:27700"` (British National Grid).
    pub crs: String,
    /// Units of the stored values. Must be `"nT"`.
    pub units: String,
    /// Height of the map surface above ground or datum (m), if known.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub altitude_m: Option<f64>,
    /// Free-text provenance, e.g. the survey name and processing steps.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
}

/// The JSON sidecar: grid geometry plus the descriptive metadata, flattened into one object.
#[derive(Serialize, Deserialize)]
struct Sidecar {
    x0: f64,
    y0: f64,
    dx: f64,
    dy: f64,
    #[serde(flatten)]
    metadata: MapMetadata,
}

/// A regular grid of magnetic anomaly values (nT) with bilinear interpolation.
///
/// The fields are private so that every `MapGrid` that exists has passed the checks in
/// [`MapGrid::new`]: positive spacing and at least 2 x 2 nodes.
#[derive(Debug, Clone, PartialEq)]
pub struct MapGrid {
    x0: f64,
    y0: f64,
    dx: f64,
    dy: f64,
    values: Array2<f64>,
}

impl MapGrid {
    /// Builds a grid from its south-west node position, spacing and `(ny, nx)` values.
    pub fn new(x0: f64, y0: f64, dx: f64, dy: f64, values: Array2<f64>) -> Result<Self, MapError> {
        if !(dx.is_finite() && dy.is_finite() && dx > 0.0 && dy > 0.0) {
            return Err(MapError::InvalidSpacing { dx, dy });
        }
        if !(x0.is_finite() && y0.is_finite()) {
            return Err(MapError::InvalidOrigin { x0, y0 });
        }
        let (ny, nx) = values.dim();
        if nx < 2 || ny < 2 {
            return Err(MapError::TooSmall { nx, ny });
        }
        Ok(Self {
            x0,
            y0,
            dx,
            dy,
            values,
        })
    }

    /// Builds a grid by sampling `f(x, y)` at every node. Handy for analytic test fields.
    ///
    /// `f` is any closure or function taking `(x, y)` and returning a value; the generic
    /// parameter `F: Fn(f64, f64) -> f64` lets the compiler inline it.
    pub fn from_fn<F>(
        x0: f64,
        y0: f64,
        dx: f64,
        dy: f64,
        nx: usize,
        ny: usize,
        f: F,
    ) -> Result<Self, MapError>
    where
        F: Fn(f64, f64) -> f64,
    {
        let values = Array2::from_shape_fn((ny, nx), |(row, col)| {
            f(x0 + col as f64 * dx, y0 + row as f64 * dy)
        });
        Self::new(x0, y0, dx, dy, values)
    }

    /// Number of columns (x direction).
    pub fn nx(&self) -> usize {
        self.values.ncols()
    }

    /// Number of rows (y direction).
    pub fn ny(&self) -> usize {
        self.values.nrows()
    }

    /// x position of the south-west node (m).
    pub fn x0(&self) -> f64 {
        self.x0
    }

    /// y position of the south-west node (m).
    pub fn y0(&self) -> f64 {
        self.y0
    }

    /// Grid spacing in x (m).
    pub fn dx(&self) -> f64 {
        self.dx
    }

    /// Grid spacing in y (m).
    pub fn dy(&self) -> f64 {
        self.dy
    }

    /// Read-only view of the stored values, shape `(ny, nx)`, row 0 = south.
    pub fn values(&self) -> &Array2<f64> {
        &self.values
    }

    /// x position of the eastern edge (last column of nodes).
    pub fn x_max(&self) -> f64 {
        self.x0 + (self.nx() - 1) as f64 * self.dx
    }

    /// y position of the northern edge (last row of nodes).
    pub fn y_max(&self) -> f64 {
        self.y0 + (self.ny() - 1) as f64 * self.dy
    }

    /// Bilinear interpolation at `(x, y)`.
    ///
    /// Returns `None` outside the grid (edges are inside) or if any of the four surrounding
    /// nodes is `NaN` (no data).
    pub fn interpolate(&self, x: f64, y: f64) -> Option<f64> {
        // Position in "grid units": fx = 2.5 means halfway between columns 2 and 3.
        let fx = (x - self.x0) / self.dx;
        let fy = (y - self.y0) / self.dy;
        let max_fx = (self.nx() - 1) as f64;
        let max_fy = (self.ny() - 1) as f64;
        // `contains` is false for NaN, so non-finite positions also return None.
        if !(0.0..=max_fx).contains(&fx) || !(0.0..=max_fy).contains(&fy) {
            return None;
        }

        // Lower-left node of the containing cell. On the far edge (fx == max_fx) we use the
        // last cell, so col + 1 stays in bounds.
        let col = (fx.floor() as usize).min(self.nx() - 2);
        let row = (fy.floor() as usize).min(self.ny() - 2);
        let tx = fx - col as f64;
        let ty = fy - row as f64;

        let v = &self.values;
        let f00 = v[[row, col]];
        let f10 = v[[row, col + 1]];
        let f01 = v[[row + 1, col]];
        let f11 = v[[row + 1, col + 1]];

        let value = (1.0 - tx) * (1.0 - ty) * f00
            + tx * (1.0 - ty) * f10
            + (1.0 - tx) * ty * f01
            + tx * ty * f11;
        // NaN in any corner propagates into `value`; report that as "no data".
        if value.is_nan() { None } else { Some(value) }
    }

    /// Loads `foo.npy` and its sidecar `foo.json` (see the module docs for the format).
    pub fn load(npy_path: impl AsRef<Path>) -> Result<(Self, MapMetadata), MapError> {
        let npy_path = npy_path.as_ref();
        let json_path = npy_path.with_extension("json");

        // `?` returns early with the error if a step fails. `map_err` first wraps the
        // library's error in our own variant so the message says which file was involved.
        let json = std::fs::read_to_string(&json_path).map_err(|source| MapError::Io {
            path: json_path.clone(),
            source,
        })?;
        let sidecar: Sidecar = serde_json::from_str(&json).map_err(|source| MapError::Json {
            path: json_path.clone(),
            source,
        })?;
        if sidecar.metadata.units != "nT" {
            return Err(MapError::WrongUnits(sidecar.metadata.units));
        }

        let values: Array2<f64> = read_npy(npy_path).map_err(|source| MapError::ReadNpy {
            path: npy_path.to_path_buf(),
            source,
        })?;
        let grid = Self::new(sidecar.x0, sidecar.y0, sidecar.dx, sidecar.dy, values)?;
        Ok((grid, sidecar.metadata))
    }

    /// Writes the grid to `foo.npy` plus the sidecar `foo.json`.
    pub fn save(&self, npy_path: impl AsRef<Path>, metadata: &MapMetadata) -> Result<(), MapError> {
        let npy_path = npy_path.as_ref();
        let json_path = npy_path.with_extension("json");

        write_npy(npy_path, &self.values).map_err(|source| MapError::WriteNpy {
            path: npy_path.to_path_buf(),
            source,
        })?;
        let sidecar = Sidecar {
            x0: self.x0,
            y0: self.y0,
            dx: self.dx,
            dy: self.dy,
            metadata: metadata.clone(),
        };
        // Serialising plain numbers and strings cannot fail, so `expect` documents that.
        let json = serde_json::to_string_pretty(&sidecar).expect("sidecar is always valid JSON");
        std::fs::write(&json_path, json).map_err(|source| MapError::Io {
            path: json_path,
            source,
        })
    }
}

impl FieldMap for MapGrid {
    fn field_at(&self, x: f64, y: f64) -> Option<f64> {
        self.interpolate(x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::{RngExt, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;

    // Origins and spacings are powers-of-two multiples so node coordinates are exact in
    // floating point. That lets the node test demand bit-exact equality.
    const X0: f64 = -1000.0;
    const Y0: f64 = 500.0;
    const DX: f64 = 250.0;
    const DY: f64 = 125.0;
    // Deliberately non-square, so swapped x/y axes cannot pass by accident.
    const NX: usize = 9;
    const NY: usize = 5;

    fn grid_of<F: Fn(f64, f64) -> f64>(f: F) -> MapGrid {
        MapGrid::from_fn(X0, Y0, DX, DY, NX, NY, f).unwrap()
    }

    /// Random points strictly inside the grid, reproducible from `seed`.
    fn interior_points(grid: &MapGrid, n: usize, seed: u64) -> Vec<(f64, f64)> {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        (0..n)
            .map(|_| {
                (
                    rng.random_range(grid.x0()..grid.x_max()),
                    rng.random_range(grid.y0()..grid.y_max()),
                )
            })
            .collect()
    }

    #[test]
    fn returns_stored_values_exactly_at_nodes() {
        let grid = grid_of(|x, y| (x * 0.003).sin() + (y * 0.001).cos());
        for row in 0..NY {
            for col in 0..NX {
                let x = X0 + col as f64 * DX;
                let y = Y0 + row as f64 * DY;
                assert_eq!(grid.interpolate(x, y), Some(grid.values()[[row, col]]));
            }
        }
    }

    #[test]
    fn linear_field_is_reproduced_exactly() {
        let f = |x: f64, y: f64| 42.0 + 0.3 * x - 0.7 * y;
        let grid = grid_of(f);
        for (x, y) in interior_points(&grid, 1000, 1) {
            assert_relative_eq!(
                grid.interpolate(x, y).unwrap(),
                f(x, y),
                max_relative = 1e-9
            );
        }
    }

    #[test]
    fn bilinear_field_is_reproduced_exactly() {
        // Bilinear interpolation is exact for a + bx + cy + dxy, including the xy cross term.
        let f = |x: f64, y: f64| 5.0 + 1e-4 * x * y;
        let grid = grid_of(f);
        for (x, y) in interior_points(&grid, 1000, 2) {
            assert_relative_eq!(
                grid.interpolate(x, y).unwrap(),
                f(x, y),
                max_relative = 1e-9
            );
        }
    }

    /// Largest interpolation error for f = x^2 at cell midpoints along x.
    fn quadratic_midpoint_error(dx: f64) -> f64 {
        let nx = ((2000.0 / dx) as usize) + 1;
        let grid = MapGrid::from_fn(0.0, 0.0, dx, dx, nx, 3, |x, _| x * x).unwrap();
        (0..nx - 1)
            .map(|col| {
                let x = (col as f64 + 0.5) * dx;
                (grid.interpolate(x, dx).unwrap() - x * x).abs()
            })
            .fold(0.0, f64::max)
    }

    #[test]
    fn quadratic_error_obeys_analytic_bound() {
        // Linear interpolation error bound: max|f''| h^2 / 8 = 2 h^2 / 8 = h^2 / 4 for f = x^2,
        // attained exactly at the midpoint of each cell.
        let grid = grid_of(|x, _| x * x);
        let bound = DX * DX / 4.0;
        for (x, y) in interior_points(&grid, 1000, 3) {
            let err = (grid.interpolate(x, y).unwrap() - x * x).abs();
            assert!(
                err <= bound * (1.0 + 1e-12),
                "error {err} exceeds bound {bound}"
            );
        }
        assert_relative_eq!(quadratic_midpoint_error(DX), bound, max_relative = 1e-9);
    }

    #[test]
    fn quadratic_error_converges_at_second_order() {
        let ratio = quadratic_midpoint_error(100.0) / quadratic_midpoint_error(50.0);
        assert!(
            (ratio - 4.0).abs() < 0.04,
            "halving dx changed error by {ratio}, want 4"
        );
    }

    #[test]
    fn outside_is_none_and_edges_are_inside() {
        let grid = grid_of(|x, y| x + y);
        assert_eq!(
            grid.interpolate(grid.x_max(), grid.y_max()),
            Some(grid.values()[[NY - 1, NX - 1]])
        );
        assert!(grid.interpolate(X0, Y0).is_some());
        assert!(grid.interpolate(X0 - 1e-6, Y0).is_none());
        assert!(grid.interpolate(X0, Y0 - 1e-6).is_none());
        assert!(grid.interpolate(grid.x_max() + 1e-6, Y0).is_none());
        assert!(grid.interpolate(X0, grid.y_max() + 1e-6).is_none());
        assert!(grid.interpolate(f64::NAN, Y0).is_none());
    }

    #[test]
    fn row_zero_is_south_and_columns_are_x() {
        // f depends on x and y with very different weights, so a swapped axis or a flipped
        // row order gives a different number.
        let grid = grid_of(|x, y| x + 1000.0 * y);
        assert_eq!(grid.values().dim(), (NY, NX));
        assert_eq!(grid.values()[[0, 0]], X0 + 1000.0 * Y0);
        assert_eq!(grid.values()[[NY - 1, 0]], X0 + 1000.0 * grid.y_max());
        assert_eq!(grid.values()[[0, NX - 1]], grid.x_max() + 1000.0 * Y0);
    }

    #[test]
    fn nan_nodes_are_reported_as_no_data() {
        // 3 rows x 4 columns, NaN at row 1, column 1.
        let mut values = Array2::from_elem((3, 4), 1.0);
        values[[1, 1]] = f64::NAN;
        let grid = MapGrid::new(0.0, 0.0, 1.0, 1.0, values).unwrap();
        // Any point whose cell has a NaN corner is "no data", even exactly on a valid node:
        // the NaN corner gets weight 0, but 0 * NaN is still NaN.
        assert!(grid.interpolate(0.5, 0.5).is_none());
        assert!(grid.interpolate(0.0, 0.0).is_none());
        // The cell spanning columns 2..3, rows 0..1 has no NaN corner.
        assert_eq!(grid.interpolate(2.5, 0.5), Some(1.0));
    }

    #[test]
    fn rejects_bad_spacing_and_tiny_grids() {
        let ok = Array2::zeros((2, 2));
        assert!(matches!(
            MapGrid::new(0.0, 0.0, -1.0, 1.0, ok.clone()),
            Err(MapError::InvalidSpacing { .. })
        ));
        assert!(matches!(
            MapGrid::new(0.0, 0.0, 1.0, f64::NAN, ok.clone()),
            Err(MapError::InvalidSpacing { .. })
        ));
        assert!(matches!(
            MapGrid::new(f64::INFINITY, 0.0, 1.0, 1.0, ok),
            Err(MapError::InvalidOrigin { .. })
        ));
        assert!(matches!(
            MapGrid::new(0.0, 0.0, 1.0, 1.0, Array2::zeros((1, 5))),
            Err(MapError::TooSmall { nx: 5, ny: 1 })
        ));
    }

    fn bng_metadata() -> MapMetadata {
        MapMetadata {
            crs: "EPSG:27700".into(),
            units: "nT".into(),
            altitude_m: Some(305.0),
            source: Some("unit test".into()),
        }
    }

    #[test]
    fn save_then_load_round_trips_exactly() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map.npy");
        let grid = grid_of(|x, y| (x * 0.003).sin() * 100.0 + y * 0.01);
        grid.save(&path, &bng_metadata()).unwrap();

        let (loaded, metadata) = MapGrid::load(&path).unwrap();
        assert_eq!(loaded, grid); // values and geometry, bit for bit
        assert_eq!(metadata, bng_metadata());
    }

    #[test]
    fn optional_metadata_fields_may_be_absent() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map.npy");
        write_npy(&path, &Array2::<f64>::zeros((3, 4))).unwrap();
        let json = r#"{"x0": 0, "y0": 0, "dx": 10, "dy": 10, "crs": "EPSG:27700", "units": "nT"}"#;
        std::fs::write(dir.path().join("map.json"), json).unwrap();

        let (grid, metadata) = MapGrid::load(&path).unwrap();
        assert_eq!((grid.nx(), grid.ny()), (4, 3));
        assert_eq!(metadata.altitude_m, None);
    }

    #[test]
    fn load_reports_bad_inputs() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("map.npy");
        let json_path = dir.path().join("map.json");

        // Missing sidecar.
        write_npy(&path, &Array2::<f64>::zeros((3, 3))).unwrap();
        assert!(matches!(MapGrid::load(&path), Err(MapError::Io { .. })));

        // Negative spacing.
        let json = r#"{"x0": 0, "y0": 0, "dx": -10, "dy": 10, "crs": "EPSG:27700", "units": "nT"}"#;
        std::fs::write(&json_path, json).unwrap();
        assert!(matches!(
            MapGrid::load(&path),
            Err(MapError::InvalidSpacing { .. })
        ));

        // Wrong units.
        let json = r#"{"x0": 0, "y0": 0, "dx": 10, "dy": 10, "crs": "EPSG:27700", "units": "T"}"#;
        std::fs::write(&json_path, json).unwrap();
        assert!(matches!(MapGrid::load(&path), Err(MapError::WrongUnits(u)) if u == "T"));

        // A 1-D array is not a map.
        let json = r#"{"x0": 0, "y0": 0, "dx": 10, "dy": 10, "crs": "EPSG:27700", "units": "nT"}"#;
        std::fs::write(&json_path, json).unwrap();
        write_npy(&path, &ndarray::Array1::<f64>::zeros(5)).unwrap();
        assert!(matches!(
            MapGrid::load(&path),
            Err(MapError::ReadNpy { .. })
        ));

        // Missing required field.
        std::fs::write(&json_path, r#"{"x0": 0, "y0": 0, "dx": 10}"#).unwrap();
        assert!(matches!(MapGrid::load(&path), Err(MapError::Json { .. })));
    }
}
