# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy", "pandas", "scipy", "verde", "harmonica"]
# ///
"""National (GB) field grids at constant altitude: equivalent sources fitted tile by tile.

The method is that of upward_bgs.py (same parameters, no retuning). Great Britain is cut into
100 km x 100 km tiles aligned to multiples of 100 km in BNG. Each tile is fitted on the data
in its core plus a 20 km halo and predicts its core nodes, which are written into one national
array with the extent and node layout of bgs_250m.npy.

Checks (pass criteria fixed before the first run):
  1. Method check per tile (tiles with >= 2000 core points): the synthetic harmonic field of
     upward_bgs.py, sampled at the tile's real positions and heights. Relative RMS <= 5% at
     every altitude. Expected to fail in most tiles, as it did over Mull; the per-tile
     numbers are the method-error table that Phase 3 uses.
  2. Hold-out per tile (fold-0 lines): ES RMS <= 1.2x cubic RMS in each tile with >= 500
     held-out points. Overall pass: >= 90% of those tiles.
  3. Seams: each tile also predicts a 2 km rim outside its core. Over the strip shared by two
     adjacent tiles, relative RMS difference <= 2% at every altitude, for every seam with
     >= 100 usable nodes.
  4. Agreement with the single-fit Mull grids (bgs_es_mull_*), 10 km border excluded:
     relative RMS difference <= 3% at every altitude.

Outputs (gitignored):

    data/grid/bgs_es_gb_<alt>m_aod.npy / .json   field at <alt> m above OD (nT)
    data/grid/bgs_es_gb_tiles.csv                per-tile checks

    uv run scripts/upward_gb.py [--altitudes 1000 2000 4000] [--tile-km 100] [--halo-km 20]
                                [--region W E S N]
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd

from bgs_common import ROOT, load_points, repair_heights
from grid_bgs import Grid, grid_region
from holdout_bgs import whole_line_key
from upward_bgs import (MULL_BOX, SOURCE, SPACING, fit, holdout_check, synthetic_errors,
                        usable)

RIM_M = 2_000.0
MIN_POINTS_CHECKS = 2_000
MIN_HOLDOUT_POINTS = 500
MIN_SEAM_NODES = 100
MULL_BORDER_M = 10_000.0


@dataclass
class Block:
    """Predictions over a tile's core plus rim: rows r0.., columns c0.. of the national grid."""

    r0: int
    c0: int
    values: dict[float, np.ndarray]  # altitude -> (rows, cols)


def index_range(axis: np.ndarray, lo: float, hi: float) -> tuple[int, int]:
    """[i0, i1) of the nodes with lo <= axis < hi."""
    return int(np.searchsorted(axis, lo, "left")), int(np.searchsorted(axis, hi, "left"))


def rel_rms(a: np.ndarray, b: np.ndarray) -> float:
    """RMS of a - b relative to the RMS of b."""
    return float(np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b**2)))


def seam_errors(p: Block, q: Block, alts: list[float]) -> tuple[int, list[float]]:
    """Relative RMS difference of two tiles' predictions where their blocks overlap."""
    def window(b: Block, alt: float, r: tuple[int, int], c: tuple[int, int]) -> np.ndarray:
        return b.values[alt][r[0] - b.r0:r[1] - b.r0, c[0] - b.c0:c[1] - b.c0]

    shape = next(iter(p.values.values())).shape
    shape_q = next(iter(q.values.values())).shape
    r = (max(p.r0, q.r0), min(p.r0 + shape[0], q.r0 + shape_q[0]))
    c = (max(p.c0, q.c0), min(p.c0 + shape[1], q.c0 + shape_q[1]))
    n, errors = 0, []
    for alt in alts:
        a, b = window(p, alt, r, c), window(q, alt, r, c)
        both = np.isfinite(a) & np.isfinite(b)
        n = int(both.sum())
        a, b = a[both], b[both]
        errors.append(float(np.sqrt(np.mean((a - b) ** 2) / np.mean(((a + b) / 2) ** 2)))
                      if n else np.nan)
    return n, errors


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--altitudes", type=float, nargs="+", default=[1000.0, 2000.0, 4000.0])
    parser.add_argument("--tile-km", type=float, default=100.0)
    parser.add_argument("--halo-km", type=float, default=20.0)
    parser.add_argument("--region", type=float, nargs=4, metavar=("W", "E", "S", "N"),
                        help="only tiles whose core overlaps this box (BNG metres; for testing)")
    args = parser.parse_args()
    alts, tile, halo = args.altitudes, args.tile_km * 1e3, args.halo_km * 1e3

    d_all = repair_heights(load_points())
    all_lines = whole_line_key(d_all).unique()
    x0, x1, y0, y1 = grid_region(d_all, SPACING)
    xs = x0 + SPACING * np.arange(round((x1 - x0) / SPACING) + 1)
    ys = y0 + SPACING * np.arange(round((y1 - y0) / SPACING) + 1)
    national = {alt: np.full((len(ys), len(xs)), np.nan) for alt in alts}
    e, n = d_all.grid_east.to_numpy(), d_all.grid_north.to_numpy()
    print(f"{len(d_all)} points, grid {len(ys)} x {len(xs)} nodes, {tile / 1e3:.0f} km tiles, "
          f"{halo / 1e3:.0f} km halo", flush=True)

    blocks: dict[tuple[int, int], Block] = {}
    rows = []
    start = time.perf_counter()
    for ty in range(int(y0 // tile), int(y1 // tile) + 1):
        for tx in range(int(x0 // tile), int(x1 // tile) + 1):
            w, s = tx * tile, ty * tile
            if args.region and not (w < args.region[1] and w + tile > args.region[0]
                                    and s < args.region[3] and s + tile > args.region[2]):
                continue
            in_core = (e >= w) & (e < w + tile) & (n >= s) & (n < s + tile)
            if not in_core.any():
                continue
            t0 = time.perf_counter()
            fit_region = (w - halo, w + tile + halo, s - halo, s + tile + halo)
            d = d_all[(e >= fit_region[0]) & (e < fit_region[1])
                      & (n >= fit_region[2]) & (n < fit_region[3])]
            label = f"tile E{w / 1e3:.0f} N{s / 1e3:.0f}"
            print(f"{label}: {int(in_core.sum())} core points, {len(d)} fitted", flush=True)

            # Grid: core plus rim, predicted from this tile's fit.
            es = fit(d)
            r = index_range(ys, s - RIM_M, s + tile + RIM_M)
            c = index_range(xs, w - RIM_M, w + tile + RIM_M)
            east, north = np.meshgrid(xs[c[0]:c[1]], ys[r[0]:r[1]])
            block = Block(r[0], c[0], {})
            row = {"tile_e_km": w / 1e3, "tile_n_km": s / 1e3, "core_points": int(in_core.sum()),
                   "fit_points": len(d), "surveys": " ".join(sorted(d_all.survey_area[in_core].unique()))}
            rc, cc = index_range(ys, s, s + tile), index_range(xs, w, w + tile)
            for alt in alts:
                ok = usable(d, east, north, alt)
                values = np.full(east.shape, np.nan)
                values[ok] = es.predict((east[ok], north[ok], np.full(ok.sum(), alt)))
                block.values[alt] = values
                national[alt][rc[0]:rc[1], cc[0]:cc[1]] = values[rc[0] - r[0]:rc[1] - r[0],
                                                                 cc[0] - c[0]:cc[1] - c[0]]
            blocks[(tx, ty)] = block

            # Checks 1 and 2.
            if in_core.sum() >= MIN_POINTS_CHECKS:
                errors = synthetic_errors(d, fit_region, alts, 500.0, "bgs")
                for alt, err in zip(alts, errors):
                    row[f"synthetic_{alt:.0f}m"] = err
                row["synthetic_pass"] = all(err <= 0.05 for err in errors)
                print("  synthetic: " + ", ".join(f"{a:.0f} m: {x:.1%}" for a, x in zip(alts, errors)))
                passed, rms = holdout_check(d, all_lines, fit_region)
                row.update(holdout_n=rms["n"], holdout_es_nt=rms["es"], holdout_cubic_nt=rms["cubic"])
                if rms["n"] >= MIN_HOLDOUT_POINTS:
                    row["holdout_pass"] = passed
            rows.append(row)
            print(f"  {time.perf_counter() - t0:.0f} s", flush=True)
    print(f"{len(blocks)} tiles in {(time.perf_counter() - start) / 60:.1f} min")

    out_dir = ROOT / "data" / "grid"
    for alt in alts:
        path = out_dir / f"bgs_es_gb_{alt:.0f}m_aod.npy"
        source = (SOURCE.format(alt=alt) + f"; fitted in {tile / 1e3:.0f} km tiles with a "
                  f"{halo / 1e3:.0f} km data halo")
        Grid(x0, y0, SPACING, national[alt]).save(path, "nT", source, altitude_m=alt)
        print(f"wrote {path.name}: {np.isfinite(national[alt]).mean():.0%} of nodes usable")
    tiles = pd.DataFrame(rows)
    for col in ["synthetic_pass", "holdout_pass", "holdout_es_nt", "holdout_cubic_nt",
                *(f"synthetic_{alt:.0f}m" for alt in alts)]:
        if col not in tiles:
            tiles[col] = np.nan
    tiles.to_csv(out_dir / "bgs_es_gb_tiles.csv", index=False, float_format="%.4g")

    print("\n1. method check per tile (criterion <= 5% at every altitude)")
    checked = tiles.dropna(subset=["synthetic_pass"])
    for alt in alts:
        col = checked[f"synthetic_{alt:.0f}m"]
        print(f"  {alt:.0f} m: median {col.median():.1%}, range {col.min():.1%}-{col.max():.1%}")
    ok1 = bool(checked.synthetic_pass.all())
    print(f"  {int(checked.synthetic_pass.sum())}/{len(checked)} tiles pass -> "
          f"{'pass' if ok1 else 'FAIL'}")

    print("2. hold-out per tile (criterion ES <= 1.2x cubic in >= 90% of tiles)")
    held = tiles.dropna(subset=["holdout_pass"])
    frac = held.holdout_pass.astype(bool).mean()
    ok2 = bool(frac >= 0.9)
    print(f"  {int(held.holdout_pass.sum())}/{len(held)} tiles pass ({frac:.0%}); median ES "
          f"{held.holdout_es_nt.median():.0f} nT, cubic {held.holdout_cubic_nt.median():.0f} nT "
          f"-> {'pass' if ok2 else 'FAIL'}")

    print("3. seams (criterion <= 2% at every altitude)")
    seams = []
    for (tx, ty), p in blocks.items():
        for q_key in [(tx + 1, ty), (tx, ty + 1)]:
            if q_key in blocks:
                count, errors = seam_errors(p, blocks[q_key], alts)
                if count >= MIN_SEAM_NODES:
                    seams.append(errors)
    seams = np.array(seams).reshape(-1, len(alts))
    ok3 = bool((seams <= 0.02).all())
    for i, alt in enumerate(alts):
        if len(seams):
            print(f"  {alt:.0f} m: median {np.median(seams[:, i]):.2%}, worst {seams[:, i].max():.2%}")
    print(f"  {int((seams <= 0.02).all(axis=1).sum())}/{len(seams)} seams pass -> "
          f"{'pass' if ok3 else 'FAIL'}")

    print("4. agreement with the single-fit Mull grids (criterion <= 3%)")
    ok4 = True
    for alt in alts:
        mull = np.load(out_dir / f"bgs_es_mull_{alt:.0f}m_aod.npy")
        mx = MULL_BOX[0] + SPACING * np.arange(mull.shape[1])
        my = MULL_BOX[2] + SPACING * np.arange(mull.shape[0])
        inner_c = (mx >= MULL_BOX[0] + MULL_BORDER_M) & (mx <= MULL_BOX[1] - MULL_BORDER_M)
        inner_r = (my >= MULL_BOX[2] + MULL_BORDER_M) & (my <= MULL_BOX[3] - MULL_BORDER_M)
        c0, r0 = int(np.searchsorted(xs, mx[0])), int(np.searchsorted(ys, my[0]))
        gb = national[alt][r0:r0 + mull.shape[0], c0:c0 + mull.shape[1]]
        a, b = gb[np.ix_(inner_r, inner_c)], mull[np.ix_(inner_r, inner_c)]
        both = np.isfinite(a) & np.isfinite(b)
        err = rel_rms(a[both], b[both])
        ok4 &= err <= 0.03
        print(f"  {alt:.0f} m: {err:.2%} over {int(both.sum())} nodes")
    print(f"  -> {'pass' if ok4 else 'FAIL'}")
    print("overall:", "PASS" if ok1 and ok2 and ok3 and ok4 else "FAIL")


if __name__ == "__main__":
    main()
