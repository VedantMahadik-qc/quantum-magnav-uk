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
| 1 | Rust workspace, map interpolation, particle filter on a synthetic map, benchmark | in progress |
| 2 | Real BGS map, sensor models, error-budget sweeps | – |
| 3 | UK navigability heatmap | – |
| 4 | Python bindings (PyO3 + maturin), plots | – |

## Layout

```
crates/magnav-core   library: maps, simulation, particle filter
crates/magnav-cli    `magnav` command-line driver
crates/magnav-py     placeholder for Python bindings (Phase 4)
docs/                research notes and citations
```

## Building and testing

```sh
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check
```

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
[BGS OGC API](https://ogcapi.bgs.ac.uk/collections/aeromag). The licence is being confirmed (see
research notes). If it is the Open Government Licence: "Contains NERC materials ©NERC [year]".
