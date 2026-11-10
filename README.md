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
| 2e | Filter estimates the INS velocity bias (marginalised particle filter), NEES consistency metric | done |
| 2e+ | Realistic INS: Schuler-loop errors from HG4930 datasheet values; filter tuned for it | done |
| 2b | Field at constant altitude via equivalent sources; Mull box done, national grids pending | partial |
| 2c, d, f | Sensor/calibration models, separate truth/filter maps, error-budget sweeps | – |
| 3 | UK navigability heatmap | – |
| 4 | Python bindings (PyO3 + maturin), plots | – |

## Layout

```
crates/magnav-core   library: maps, simulation, particle filter
crates/magnav-cli    `magnav` command-line driver
crates/magnav-py     placeholder for Python bindings (Phase 4)
scripts/             data pipeline (Python, run with uv): fetch, report, grid, validate, upward-continue
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
uv run scripts/upward_bgs.py           # Mull box at 1, 2, 4 km above sea level (~1 min)
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
case). Each run is 60 km east, and the INS ends 1000 m off because of a 2 m/s velocity bias.
Lines are chosen by `scripts/pick_demo_lines.py`. Final errors, 3 seeds:

| Line | RMS gradient | Position-only filter | With bias estimation |
|---|---|---|---|
| Isle of Mull (E 144 km, N 730 km) | 893 nT/km | 36–37 m (predicted < 200 m) | 1.8–2.3 m (predicted < 10 m) |
| Dorset chalk (E 324 km, N 100 km) | 1.4 nT/km | 677–680 m (predicted ≥ 800 m: **wrong**) | 26–69 m, error 0.14–0.37× reported spread (predicted ≤ 3×) |

The first Dorset prediction failed in an instructive way. The position-only filter got a fix
early (76 m at 100 s), then drifted with the INS bias while reporting only ~110 m of spread,
ending about 6σ off: *overconfident*. Estimating the bias fixes both problems, and the filter
learns the 2 m/s bias even over flat chalk.

How the bias is estimated matters. Simply adding it to each particle's state failed the
pre-set bias check: resampling collapses the bias hypotheses, and the filter ends confidently
wrong (sample impoverishment). The filter now *marginalises* the bias instead. Each particle carries an exact Kalman
estimate of it (Schön, Gustafsson & Nordlund 2005), which is also the method used in the
MagNav literature (Canciani & Raquet 2016, as cited in Sengupta 2025).

These runs use a perfect sensor, so the filter is *underconfident*. NEES, the normalised
estimation error squared, is about 0.01–0.1 against an ideal of 2. I first blamed the 5 m/√s
position process noise. A sweep showed it isn't the cause: the filter's 2 nT likelihood meets
readings that are actually exact. See the next section.

## A realistic INS (HG4930)

The 2 m/s constant bias above is a toy. An unaided tactical-grade INS is much worse, and its
error grows with time. `crates/magnav-core/src/ins.rs` models each horizontal axis with the
standard Schuler loop (velocity error and tilt coupled through gravity and Earth radius),
propagated exactly. It uses the **Honeywell HG4930** (grade CA51) values from Honeywell's
*Performance and Environmental Information* (N61-1775-000-000, 06/17): gyro bias 1.0 °/h,
ARW 0.09 °/√h, accelerometer bias 0.3 mg, VRW 0.09 m/s/√h.
After a 10-minute GNSS outage this INS is typically **0.2–4 km** off, mostly from the gyro
bias tilting the platform.

The filter's bias model is sized from the same IMU model. With 2 nT sensor noise matching the
filter's likelihood, it is consistent: mean NEES is 1.2–1.3 on synthetic maps for every
q in 0.5–5 m/√s, and q = 1 m/√s is used. On the real map, 60 km east, starting on the true
position at GNSS loss, 3 seeds:

| Line | INS alone | MagNav filter | Error ÷ reported spread |
|---|---|---|---|
| Isle of Mull | 214–3949 m | **4.5–6.6 m** | 0.26–0.37 |
| Dorset chalk | 214–3949 m | **283–1132 m** | 0.98–2.98 |

So where the crust is strongly magnetic, MagNav holds a realistic INS to metres. Over the
chalk it only roughly halves the error. On one seed it is slightly worse than an INS that
happened to draw small biases. This is the navigability contrast Phase 3 will map
nationally. (Predictions fixed before running: Mull < 50 m, and error ≤ 3× spread on both lines.
Both passed; Dorset seed 3 only just, at 2.98×.)

```sh
cargo run --release -p magnav-cli -- demo --map-file data/grid/bgs_250m.npy \
    --start-x 144000 --start-y 730000 --ins hg4930 --sensor-noise-nt 2
```

## Altitude (Phase 2b)

The mixed-height survey is levelled and continued upward with equivalent sources
([docs/upward-continuation.md](docs/upward-continuation.md)); so far only the Mull box is done. Over Mull, a realistic HG4930 INS
(median of 3 seeds, 60 km flight):

| Flying at | MagNav final error |
|---|---|
| 305 m above ground (draped) | 5.8 m |
| 1000 m above sea level | 16.5 m |
| 2000 m above sea level | 74 m |
| 4000 m above sea level | 310 m |

The prediction that error grows with altitude at every step held. Altitude costs about 50×
here, which matches the sharp drop in field gradient: 87 → 16 nT/km from 1 to 4 km.

The method check failed in an informative way. Equivalent sources recover a known synthetic
field to 2–3% when it is smooth or sampled without gaps. Sampled like the real survey, the
error at 1000 m is 13%: the 2 km line spacing never measured the short wavelengths between
lines. That is a hard limit of the 1960s data at low altitude, not of the processing.
At survey level, equivalent sources still predict held-out lines better than the cubic grid
(119 against 145 nT RMS).

## Data findings

- Three surveys in southern and eastern England were flown at a constant barometric height
  (549 m or 457 m above sea level). All others, including all of Scotland, followed the
  terrain at 152 m or 305 m. So the map is not at one altitude ([docs/data-report.md](docs/data-report.md)).
- At a 500 m grid spacing the map failed its own fidelity check (18 nT RMS against ≤ 10 nT):
  500 m blocks average away the steep Scottish anomalies. 250 m, the median spacing of points
  along a line, passes narrowly (9.5 nT).

## What is validated so far (synthetic maps)

Each check has its pass criterion fixed before it was run; see the test files for details.

| Check | Criterion | Result |
|---|---|---|
| Bilinear interpolation, linear and bilinear fields | exact to 1e-9 relative | pass |
| Bilinear interpolation, x^2 | error <= h^2/4, second-order convergence | pass |
| Filter on a strongly varying synthetic map, 5 seeds | final error < 100 m and < 0.2x INS | 23-52 m vs 1000 m INS |
| Filter on a flat map, 5 seeds | no resampling; spread = sqrt(s0^2 + q^2 t) +/-5%; tracks INS; final error >= 0.8x INS | 982-1014 m (no convergence, as required) |
| Serial vs rayon | bit-identical results on 1, 3, 16 threads | pass |
| Bias estimation on, strongly varying map, 5 seeds | final error < 0.5x position-only filter; bias within 0.2 m/s | 0.3–1.2 m (vs 23–52 m); bias within 0.045 m/s |
| Bias estimation on, consistency | mean NEES over last 300 s <= 6 | 0.01 (conservative) |
| Bias estimation on, flat map | final error >= 0.8x INS | 972–1026 m (no convergence, as required) |
| Schuler INS, bias-only responses | match closed-form solutions to 1e-3 | pass |
| Schuler INS, 2000-flight Monte Carlo | velocity-error std within 5% of the sizing formula | pass |
| HG4930 INS, q sweep, 2 nT matched noise* | best NEES in [0.5, 8]; final < 0.2x INS | q = 1: NEES 1.28; 4–32 m vs 214–3949 m |

\* Criterion revised after a diagnostic: with an ideal sensor it failed (NEES 0.17–0.19).
It is a regression check, not a prediction.

One finding from building these tests: where the map is locally flat, the filter correctly
falls back to INS dead reckoning and its error grows at the INS drift rate. That is the
mechanism behind the "where in the UK does MagNav work" question.

## Particle filter benchmark

`cargo bench -p magnav-core`. One full filter step (predict, weigh, normalise, resample when
ESS < N/2) on a 100 m synthetic grid. Measured on an AMD Ryzen 9 7940HS (8 cores / 16
threads), Windows 11, rustc 1.98.1, with criterion medians:

| Particles | Serial | Rayon (16 threads) | Speed-up |
|---:|---:|---:|---:|
| 1 000 | 46 µs (22 M particles/s) | 35 µs (29 M/s) | 1.3× |
| 10 000 | 467 µs (21 M/s) | 187 µs (54 M/s) | 2.5× |
| 100 000 | 4.61 ms (22 M/s) | 1.11 ms (90 M/s) | 4.1× |

Re-measured 2026-11-10, after adding bias estimation. The per-particle bias update runs either
way, so the numbers hold with it on or off. Serial times matched the earlier run to within
3%. The rayon times vary by about 15% between sessions on this laptop.

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
