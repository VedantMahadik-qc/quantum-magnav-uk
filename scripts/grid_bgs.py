# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "pandas", "scipy", "verde"]
# ///
"""Grid the BGS aeromag points onto a regular British National Grid raster.

Method (deliberately simple; Phase 2b compares it with equivalent sources):
1. median of the points in each `spacing` x `spacing` block (removes along-line clustering);
2. Clough-Tocher cubic interpolation on the Delaunay triangulation (verde.Cubic);
3. cells further than `--mask-km` from any raw point are set to NaN (no data).

Outputs (gitignored), all in the layout `MapGrid::load` expects (row 0 = south):

    data/grid/bgs_<spacing>m.npy / .json      anomaly (nT)
    data/grid/bgs_<spacing>m_aod.npy / .json  observation height above OD of the nearest point (m)
    data/grid/check_points.csv                random off-node points + numpy bilinear value
    data/grid/raw_sample.csv                  random raw points (x, y, nT) for the Rust test

    uv run scripts/grid_bgs.py [--spacing 250] [--mask-km 2]
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import verde as vd
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from bgs_common import ROOT, load_points

SOURCE = (
    "BGS Great Britain Aeromagnetic Survey (1955-65), mag_igrf90; block median + verde.Cubic. "
    "Mixed observation surfaces: draped 305 m (most), draped 152 m (SW England), "
    "constant 549 m AOD (CA55_SOUTH, CA57), constant 457 m AOD (HG56)."
)


@dataclass
class Grid:
    """A gridded field: `values[row, col]` sits at (x0 + col*spacing, y0 + row*spacing)."""

    x0: float
    y0: float
    spacing: float
    values: np.ndarray

    def node_axes(self) -> tuple[np.ndarray, np.ndarray]:
        ny, nx = self.values.shape
        return self.x0 + self.spacing * np.arange(nx), self.y0 + self.spacing * np.arange(ny)

    def bilinear(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Same rule as Rust `MapGrid::interpolate`: NaN if any corner is NaN or off-grid."""
        xs, ys = self.node_axes()
        interp = RegularGridInterpolator((ys, xs), self.values, bounds_error=False,
                                         fill_value=np.nan)
        return interp(np.column_stack([y, x]))

    def save(self, npy_path: Path, units: str, source: str, altitude_m: float | None) -> None:
        npy_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(npy_path, self.values.astype(np.float64))
        meta = {"x0": self.x0, "y0": self.y0, "dx": self.spacing, "dy": self.spacing,
                "crs": "EPSG:27700", "units": units, "source": source}
        if altitude_m is not None:
            meta["altitude_m"] = altitude_m
        npy_path.with_suffix(".json").write_text(json.dumps(meta, indent=2) + "\n")


def grid_region(d: pd.DataFrame, spacing: float) -> tuple[float, float, float, float]:
    """Data bounding box snapped outwards to multiples of the spacing."""
    def lo(v): return float(np.floor(v.min() / spacing) * spacing)
    def hi(v): return float(np.ceil(v.max() / spacing) * spacing)
    return lo(d.grid_east), hi(d.grid_east), lo(d.grid_north), hi(d.grid_north)


def make_grid(d: pd.DataFrame, spacing: float, mask_m: float,
              region: tuple[float, float, float, float] | None = None) -> Grid:
    """Grid the anomaly in `d`. Shared with holdout_bgs.py so both use the same method."""
    region = region or grid_region(d, spacing)
    coords = (d.grid_east.to_numpy(float), d.grid_north.to_numpy(float))
    block_coords, block_values = vd.BlockReduce("median", spacing=spacing).filter(
        coords, d.mag_igrf90.to_numpy(float))
    gridder = vd.Cubic().fit(block_coords, block_values)

    # verde's grid_coordinates returns (easting, northing) arrays of shape (ny, nx) with
    # northing increasing along axis 0: row 0 = south, which is exactly the Rust convention.
    east, north = vd.grid_coordinates(region, spacing=spacing)
    values = gridder.predict((east, north))
    far = ~vd.distance_mask(coords, maxdist=mask_m, coordinates=(east, north))
    values[far] = np.nan
    return Grid(x0=region[0], y0=region[2], spacing=spacing, values=values)


def nearest_layer(d: pd.DataFrame, column: str, like: Grid) -> Grid:
    """Value of `column` at the nearest raw point to each node (NaN where `like` is NaN)."""
    xs, ys = like.node_axes()
    east, north = np.meshgrid(xs, ys)
    tree = cKDTree(np.column_stack([d.grid_east, d.grid_north]))
    _, idx = tree.query(np.column_stack([east.ravel(), north.ravel()]))
    values = d[column].to_numpy(float)[idx].reshape(east.shape)
    values[np.isnan(like.values)] = np.nan
    return Grid(like.x0, like.y0, like.spacing, values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--spacing", type=float, default=250.0, help="grid spacing (m)")
    parser.add_argument("--mask-km", type=float, default=2.0,
                        help="NaN beyond this distance from data (km)")
    args = parser.parse_args()
    out_dir = ROOT / "data" / "grid"
    stem = f"bgs_{args.spacing:.0f}m"

    d = load_points()
    grid = make_grid(d, args.spacing, args.mask_km * 1e3)
    grid.save(out_dir / f"{stem}.npy", "nT", SOURCE, altitude_m=305.0)
    nearest_layer(d, "aod", grid).save(
        out_dir / f"{stem}_aod.npy", "m",
        "observation height above Ordnance Datum of the nearest BGS point", altitude_m=None)

    rng = np.random.default_rng(0)
    xs, ys = grid.node_axes()
    # Random positions (not nodes), so the Rust test exercises the interpolation weights.
    cx = rng.uniform(xs[0], xs[-1], 20_000)
    cy = rng.uniform(ys[0], ys[-1], 20_000)
    cv = grid.bilinear(cx, cy)
    keep = ~np.isnan(cv)
    pd.DataFrame({"x": cx[keep], "y": cy[keep], "value": cv[keep]}).head(1000).to_csv(
        out_dir / "check_points.csv", index=False, float_format="%.17g")
    d.sample(10_000, random_state=0)[["grid_east", "grid_north", "mag_igrf90"]].set_axis(
        ["x", "y", "value"], axis=1).to_csv(out_dir / "raw_sample.csv", index=False)

    # Pass criterion from the plan: training-point RMS <= 10 nT (one contour interval).
    resid = grid.bilinear(d.grid_east.to_numpy(float), d.grid_north.to_numpy(float)) - d.mag_igrf90
    resid = resid[~np.isnan(resid)]
    rms = float(np.sqrt(np.mean(resid**2)))
    ny, nx = grid.values.shape
    print(f"{stem}: {nx} x {ny} nodes, {np.isnan(grid.values).mean():.1%} NaN")
    print(f"training residual at {resid.size} raw points: RMS {rms:.1f} nT, "
          f"median |r| {np.median(np.abs(resid)):.1f} nT, "
          f"95% |r| {np.quantile(np.abs(resid), 0.95):.1f} nT")
    print("PASS" if rms <= 10.0 else "FAIL", "(criterion: training RMS <= 10 nT)")


if __name__ == "__main__":
    main()
