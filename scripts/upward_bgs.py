# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "pandas", "scipy", "verde", "harmonica"]
# ///
"""Field at constant altitude from the mixed-height BGS survey, via equivalent sources.

The survey was flown on several surfaces (draped 152/305 m, constant 457/549 m above OD; see
docs/data-report.md). Equivalent sources (harmonica.EquivalentSourcesGB) fit a layer of point
sources to every data point at its real height (`aod`), and then predict the field at any
altitude. That levels the mixed surfaces and continues the field upward in one step.

Steps (pass criteria fixed before the first run):
  1. Method check on a known field: a harmonic synthetic field (buried point masses) sampled
     at the real BGS point positions and heights. Predictions at each altitude must have
     relative RMS error <= 5%.
  2. Real data, at survey level: hold out the fold-0 lines of holdout_bgs.py and predict
     them. RMS must be <= 1.2x the cubic baseline's on the same points.
  3. Write grids at the requested altitudes (m above OD), masked to within 2 km of data and
     >= 150 m above the ground.

Parameters (depth 1000 m, damping 1, 500 m source blocks, 250 m data blocks) were chosen on a
separate set of tuning lines (seed 100) that excludes the step-2 lines.

    uv run scripts/upward_bgs.py [--region W E S N] [--altitudes 1000 2000 4000]
"""

from __future__ import annotations

import argparse

import harmonica as hm
import numpy as np
import pandas as pd
import verde as vd

from bgs_common import ROOT, load_points
from grid_bgs import Grid, make_grid, nearest_layer
from holdout_bgs import HOLDOUT_FRACTION, whole_line_key

# Isle of Mull and surroundings: contains the demo flight line (E 144-204 km at N 730 km).
MULL_BOX = (100e3, 250e3, 680e3, 780e3)
SPACING = 250.0
MASK_DISTANCE_M = 2_000.0
MIN_CLEARANCE_M = 150.0
ES_PARAMS = dict(depth=1000.0, damping=1.0, block_size=500.0, window_size=20e3, random_state=0)
SOURCE = ("BGS aeromag (1955-65), equivalent sources (harmonica EquivalentSourcesGB, depth "
          "1000 m, damping 1, 500 m blocks) fitted at each point's observation height and "
          "predicted at {alt:.0f} m above Ordnance Datum")


def xyz(d: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (d.grid_east.to_numpy(float), d.grid_north.to_numpy(float), d.aod.to_numpy(float))


def fit(d: pd.DataFrame, values: np.ndarray | None = None) -> hm.EquivalentSourcesGB:
    """Fit equivalent sources to points `d` (default values: the anomaly column)."""
    values = d.mag_igrf90.to_numpy(float) if values is None else values
    coords, v = vd.BlockReduce("median", spacing=SPACING, drop_coords=False).filter(xyz(d), values)
    return hm.EquivalentSourcesGB(**ES_PARAMS).fit(coords, v)


def altitude_mask(d: pd.DataFrame, region, altitude: float) -> tuple[tuple, np.ndarray]:
    """Grid nodes at `altitude` and a boolean mask of the usable ones."""
    east, north = vd.grid_coordinates(region, spacing=SPACING)
    near = vd.distance_mask(xyz(d)[:2], maxdist=MASK_DISTANCE_M, coordinates=(east, north))
    template = Grid(region[0], region[2], SPACING, np.where(near, 0.0, np.nan))
    ground = nearest_layer(d, "ground_elev", template).values
    ok = near & (altitude >= ground + MIN_CLEARANCE_M)
    return (east, north, np.full_like(east, altitude)), ok


def synthetic_errors(d: pd.DataFrame, region, altitudes: list[float], min_depth: float,
                     sampling: str) -> list[float]:
    """Relative RMS error at each altitude when recovering a known harmonic field.

    Point masses lie between 5 km and `min_depth` below Ordnance Datum (below all terrain).
    `sampling` is "bgs" (the real point positions and heights) or "uniform" (as many points,
    spread uniformly, each at the height of its nearest real point).
    """
    rng = np.random.default_rng(42)
    n = 300
    masses_xyz = (rng.uniform(region[0], region[1], n), rng.uniform(region[2], region[3], n),
                  rng.uniform(-5_000.0, -min_depth, n))
    masses = rng.normal(0.0, 1e10, n)

    def field(coords):
        return hm.point_gravity(coords, masses_xyz, masses, field="g_z")

    obs = d
    if sampling == "uniform":
        obs = d.copy()
        obs["grid_east"] = rng.uniform(region[0], region[1], len(d))
        obs["grid_north"] = rng.uniform(region[2], region[3], len(d))
        obs["aod"] = nearest_layer_at(d, "aod", obs.grid_east, obs.grid_north)
    es = fit(obs, field(xyz(obs)))
    errors = []
    for alt in altitudes:
        coords, ok = altitude_mask(d, region, alt)
        exact, pred = field(coords)[ok], es.predict(coords)[ok]
        errors.append(float(np.sqrt(np.mean((pred - exact) ** 2)) / np.sqrt(np.mean(exact**2))))
    return errors


def nearest_layer_at(d: pd.DataFrame, column: str, x, y) -> np.ndarray:
    from scipy.spatial import cKDTree

    _, idx = cKDTree(np.column_stack([d.grid_east, d.grid_north])).query(np.column_stack([x, y]))
    return d[column].to_numpy(float)[idx]


def synthetic_check(d: pd.DataFrame, region, altitudes: list[float]) -> bool:
    """Criterion 1 (real sampling, sources down to 0.5 km), plus two diagnostic variants."""
    ok_all = True
    for label, min_depth, sampling, is_criterion in [
        ("real sampling, sources 0.5-5 km deep", 500.0, "bgs", True),
        ("diagnostic: real sampling, sources 2.5-5 km deep", 2_500.0, "bgs", False),
        ("diagnostic: gap-free sampling, sources 0.5-5 km deep", 500.0, "uniform", False),
    ]:
        errors = synthetic_errors(d, region, altitudes, min_depth, sampling)
        cells = ", ".join(f"{a:.0f} m: {e:.1%}" for a, e in zip(altitudes, errors))
        verdict = ""
        if is_criterion:
            passed = all(e <= 0.05 for e in errors)
            ok_all &= passed
            verdict = " -> pass" if passed else " -> FAIL (criterion <= 5% at every altitude)"
        print(f"  {label}: {cells}{verdict}")
    return ok_all


def holdout_check(d: pd.DataFrame, all_lines: np.ndarray, region) -> tuple[bool, dict]:
    """Criterion 2: ES vs the cubic baseline on held-out lines, at survey level."""
    test_lines = np.random.default_rng(0).choice(
        all_lines, size=int(HOLDOUT_FRACTION * len(all_lines)), replace=False)
    is_test = whole_line_key(d).isin(test_lines).to_numpy()
    train, test = d[~is_test], d[is_test]
    truth = test.mag_igrf90.to_numpy(float)

    cubic = make_grid(train, SPACING, MASK_DISTANCE_M, region=region)
    r_cubic = cubic.bilinear(test.grid_east.to_numpy(float), test.grid_north.to_numpy(float)) - truth
    r_es = fit(train).predict(xyz(test)) - truth
    common = ~np.isnan(r_cubic)
    rms = {"cubic": float(np.sqrt(np.mean(r_cubic[common] ** 2))),
           "es": float(np.sqrt(np.mean(r_es[common] ** 2))), "n": int(common.sum())}
    passed = rms["es"] <= 1.2 * rms["cubic"]
    print(f"  hold-out ({rms['n']} points on {len(set(test_lines) & set(whole_line_key(d)))} "
          f"lines): ES RMS {rms['es']:.1f} nT, cubic RMS {rms['cubic']:.1f} nT "
          f"-> {'pass' if passed else 'FAIL'} (criterion ES <= 1.2x cubic)")
    return passed, rms


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--region", type=float, nargs=4, default=MULL_BOX,
                        metavar=("W", "E", "S", "N"), help="BNG metres")
    parser.add_argument("--altitudes", type=float, nargs="+", default=[1000.0, 2000.0, 4000.0])
    args = parser.parse_args()
    region = tuple(args.region)

    all_points = load_points()
    all_lines = whole_line_key(all_points).unique()
    d = all_points[all_points.grid_east.between(region[0], region[1])
                   & all_points.grid_north.between(region[2], region[3])]
    print(f"{len(d)} points in region {region}")

    print("1. method check on a synthetic field")
    ok1 = synthetic_check(d, region, args.altitudes)
    print("2. survey-level hold-out against the cubic grid")
    ok2, _ = holdout_check(d, all_lines, region)

    print("3. grids at altitude")
    es = fit(d)
    out_dir = ROOT / "data" / "grid"
    tag = "mull" if region == MULL_BOX else "region"
    for alt in args.altitudes:
        coords, ok = altitude_mask(d, region, alt)
        values = np.where(ok, es.predict(coords), np.nan)
        grid = Grid(region[0], region[2], SPACING, values)
        path = out_dir / f"bgs_es_{tag}_{alt:.0f}m_aod.npy"
        grid.save(path, "nT", SOURCE.format(alt=alt), altitude_m=alt)
        gy, gx = np.gradient(values, SPACING)
        grad = np.sqrt(np.nanmean(gx**2 + gy**2)) * 1e3
        print(f"  {path.name}: {ok.mean():.0%} usable, field std {np.nanstd(values):.0f} nT, "
              f"RMS gradient {grad:.0f} nT/km")
    print("overall:", "PASS" if ok1 and ok2 else "FAIL")


if __name__ == "__main__":
    main()
