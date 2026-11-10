# quantum-magnav-uk

**Where in the UK can you navigate without GPS?**

Magnetic anomaly navigation (MagNav) estimates position by matching magnetometer readings against
a map of the Earth's crustal magnetic field. This project simulates MagNav over Great Britain to
answer two questions:

1. Which error source limits position accuracy: sensor noise, sensor drift, platform calibration,
   map resolution/error, or altitude?
2. Where in the UK does magnetic navigation work at all?

Hypothesis: sensor sensitivity is **not** the bottleneck. Calibration error, map quality (the UK
national aeromagnetic survey dates from 1955–65) and altitude dominate. This is a hypothesis to
test, not a result.

## Status

| Phase | Content | State |
|---|---|---|
| 0 | Literature, data and sensor-parameter research ([docs/research-notes.md](docs/research-notes.md)) | done |
| 1 | Rust workspace, map interpolation, particle filter on a synthetic map, benchmark | done |
| 2a | Real BGS map: download, grid, hold-out validation, first real-map flights | done |
| 2b–f | Upward continuation, sensor/calibration models, velocity-bias state, error-budget sweeps | – |
| 3 | UK navigability heatmap | – |
| 4 | Python bindings (PyO3 + maturin), plots | – |

## Layout

```
crates/magnav-core   library: maps, simulation, particle filter
crates/magnav-cli    `magnav` command-line driver
crates/magnav-py     placeholder for Python bindings (Phase 4)
scripts/             data pipeline (Python, run with uv): fetch, report, grid, validate
docs/                research notes, data report, map validation
```

## Building and testing

```sh
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
cargo run --release -p magnav-cli -- demo --map gaussian   # filter vs INS on a synthetic map
cargo run --release -p magnav-cli -- demo --map flat       # same flight, no map information
```

## Getting the map data

The BGS points are downloaded, not committed. The scripts use [uv](https://docs.astral.sh/uv/)
inline dependencies, so no virtual environment needs setting up.

```sh
uv run scripts/fetch_bgs_aeromag.py    # 541,508 points -> data/raw/ (~2 min)
uv run scripts/bgs_data_report.py      # docs/data-report.md
uv run scripts/grid_bgs.py             # 250 m grid -> data/grid/bgs_250m.npy (+ .json)
uv run scripts/holdout_bgs.py          # docs/map-validation.md (~2 min)
cargo test -p magnav-core --test bgs_grid -- --ignored   # Rust reads the grid correctly
cargo run --release -p magnav-cli -- demo --map-file data/grid/bgs_250m.npy --start-x 144000 --start-y 730000
```

## First results on the real map (Phase 2a)

**How good is the 1955–65 map?** Removing 10% of flight lines and predicting them from the
rest ([docs/map-validation.md](docs/map-validation.md)) gives the error between survey lines.
This is an overestimate, because removing a line doubles the local line spacing:

| Region | Held-out RMS error |
|---|---|
| Southern and eastern England (CA57, HG56, CA55_SOUTH) | 5–9 nT |
| Rest of England, Wales and southern Scotland (CA55_NORTH, CA58–CA60, HG57, HG58, HG61) | 15–37 nT |
| Central and northern Scotland (CA63, HG62, HG64, HG65) | 83–127 nT |

For comparison, a sensor's noise floor is pT to nT, and good aircraft calibration reaches
~5 nT (25 nT conservative; Sengupta 2025). In central and northern Scotland the map error is several times larger
than either. The map is least accurate exactly where the field is strongest, i.e. where it
carries the most position information.

**Flying over it** (perfect sensor and the same map used for truth and filter, so a best
case). Each run is 60 km east, and the INS ends 1000 m off. Lines are chosen by
`scripts/pick_demo_lines.py`:

| Line | RMS gradient | Prediction (fixed before running) | Result, 3 seeds |
|---|---|---|---|
| Isle of Mull (E 144 km, N 730 km) | 893 nT/km | filter < 200 m | 36–37 m: as predicted |
| Dorset chalk (E 324 km, N 100 km) | 1.4 nT/km | filter ≥ 800 m (no fix) | 677–680 m: **prediction wrong** |

Over Dorset the filter does get a fix early (76 m at 100 s). It then drifts with the INS
velocity bias while reporting only ~110 m of spread, so it ends about 6σ off and
overconfident. The cause is that the filter state has no velocity bias. Adding it, and checking
the filter's consistency (error versus reported spread), is the next filter change.

## Data findings

- Three surveys in southern and eastern England were flown at a constant barometric height
  (549 m or 457 m above sea level). All others, including all of Scotland, followed the
  terrain at 152 m or 305 m. So the map is not at one altitude ([docs/data-report.md](docs/data-report.md)).
- At a 500 m grid spacing the map failed its own fidelity check (18 nT RMS against ≤ 10 nT):
  500 m blocks average away the steep Scottish anomalies. 250 m, the median spacing of points
  along a line, passes narrowly (9.5 nT).

## What is validated so far (Phase 1)

Each check has its pass criterion fixed before it was run; see the test files for details.

| Check | Criterion | Result |
|---|---|---|
| Bilinear interpolation, linear and bilinear fields | exact to 1e-9 relative | pass |
| Bilinear interpolation, x^2 | error <= h^2/4, second-order convergence | pass |
| Filter on a strongly varying synthetic map, 5 seeds | final error < 100 m and < 0.2x INS | 23-52 m vs 1000 m INS |
| Filter on a flat map, 5 seeds | no resampling; spread = sqrt(s0^2 + q^2 t) +/-5%; tracks INS; final error >= 0.8x INS | 982-1014 m (no convergence, as required) |
| Serial vs rayon | bit-identical results on 1, 3, 16 threads | pass |

One finding from building these tests: where the map is locally flat, the filter correctly
falls back to INS dead reckoning and its error grows at the INS drift rate. That is the
mechanism behind the "where in the UK does MagNav work" question.

## Particle filter benchmark

`cargo bench -p magnav-core`. One full filter step (predict, weigh, normalise, resample when
ESS < N/2) on a 100 m synthetic grid. Measured on an AMD Ryzen 9 7940HS (8 cores / 16
threads), Windows 11, rustc 1.98.1, with criterion medians:

| Particles | Serial | Rayon (16 threads) | Speed-up |
|---:|---:|---:|---:|
| 1 000 | 45 µs (22 M particles/s) | 38 µs (26 M/s) | 1.2× |
| 10 000 | 455 µs (22 M/s) | 231 µs (43 M/s) | 2.0× |
| 100 000 | 4.72 ms (21 M/s) | 1.29 ms (78 M/s) | 3.7× |

The speed-up stops well short of 8× because only the per-particle predict/weigh step runs
in parallel (Amdahl's law). Normalisation, the weighted sums and resampling stay serial so
that results are bit-identical to the serial filter (tested), and at small N thread overhead
dominates.

## Conventions

Positions are planar metres (x = east, y = north; British National Grid for real data), time in
seconds, magnetic field in nanotesla (nT).

## Data and licences

Map data is **not** committed. The BGS Great Britain Aeromagnetic Survey is fetched from the
[BGS OGC API](https://ogcapi.bgs.ac.uk/collections/aeromag). The collection metadata states no
licence; BGS were asked to confirm on 2026-11-10. Until then only code and summary statistics
are committed. If it is the Open Government Licence: "Contains NERC materials ©NERC [year]".
