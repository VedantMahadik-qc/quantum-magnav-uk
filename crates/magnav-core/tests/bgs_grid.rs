//! Round-trip checks of the real BGS grid written by `scripts/grid_bgs.py`.
//!
//! These need downloaded data, so they are `#[ignore]`d and CI skips them. Run with
//!
//! ```sh
//! cargo test -p magnav-core --test bgs_grid -- --ignored
//! ```
//!
//! The grid path defaults to `data/grid/bgs_250m.npy` in the workspace and can be overridden
//! with the `MAGNAV_BGS_GRID` environment variable. Pass criteria were fixed before running.

use std::fs;
use std::path::PathBuf;

use magnav_core::MapGrid;

fn grid_path() -> PathBuf {
    match std::env::var_os("MAGNAV_BGS_GRID") {
        Some(p) => PathBuf::from(p),
        // `env!` is expanded at compile time to this crate's directory, so the default works
        // whatever directory `cargo test` is started from.
        None => PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/grid/bgs_250m.npy"),
    }
}

/// Reads an `x,y,value` CSV (with header) written next to the grid.
fn read_xyz(name: &str) -> Vec<(f64, f64, f64)> {
    let path = grid_path().with_file_name(name);
    let text =
        fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    text.lines()
        .skip(1)
        .map(|line| {
            // Parse three comma-separated numbers; `collect` into a Vec, then destructure.
            let v: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
            (v[0], v[1], v[2])
        })
        .collect()
}

fn load() -> MapGrid {
    let (grid, meta) = MapGrid::load(grid_path()).expect("load BGS grid");
    assert_eq!(meta.crs, "EPSG:27700");
    grid
}

#[test]
#[ignore = "needs data/grid from scripts/grid_bgs.py"]
fn grid_lies_inside_british_national_grid() {
    let b = load().bounds();
    assert!(
        b.x_min >= 0.0 && b.x_max <= 700_000.0,
        "easting range {b:?}"
    );
    assert!(
        b.y_min >= 0.0 && b.y_max <= 1_300_000.0,
        "northing range {b:?}"
    );
}

/// Rust and numpy must implement the same bilinear rule on the same array.
#[test]
#[ignore = "needs data/grid from scripts/grid_bgs.py"]
fn rust_bilinear_matches_numpy() {
    let grid = load();
    let points = read_xyz("check_points.csv");
    assert!(
        points.len() >= 500,
        "too few check points: {}",
        points.len()
    );
    for (x, y, expected) in points {
        let got = grid
            .interpolate(x, y)
            .unwrap_or_else(|| panic!("no value at ({x}, {y}); numpy had {expected}"));
        let tol = 1e-9 * expected.abs().max(1.0);
        assert!(
            (got - expected).abs() <= tol,
            "({x}, {y}): rust {got}, numpy {expected}"
        );
    }
}

/// The grid must reproduce the raw survey points it was built from. A transposed or
/// upside-down array would give errors of hundreds of nT, so this checks orientation.
#[test]
#[ignore = "needs data/grid from scripts/grid_bgs.py"]
fn grid_reproduces_raw_points() {
    let grid = load();
    // `filter_map` keeps only the points where `interpolate` returned `Some`.
    let mut abs_err: Vec<f64> = read_xyz("raw_sample.csv")
        .into_iter()
        .filter_map(|(x, y, v)| grid.interpolate(x, y).map(|g| (g - v).abs()))
        .collect();
    assert!(
        abs_err.len() >= 9_000,
        "only {} raw points on the grid",
        abs_err.len()
    );
    // f64 is not `Ord` (NaN breaks total ordering), so sort with `total_cmp`.
    abs_err.sort_by(f64::total_cmp);
    let median = abs_err[abs_err.len() / 2];
    assert!(median <= 10.0, "median |grid - raw| = {median} nT");
}
