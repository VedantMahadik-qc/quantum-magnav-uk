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

## Conventions

Positions are planar metres (x = east, y = north; British National Grid for real data), time in
seconds, magnetic field in nanotesla (nT).

## Data and licences

Map data is **not** committed. The BGS Great Britain Aeromagnetic Survey is fetched from the
[BGS OGC API](https://ogcapi.bgs.ac.uk/collections/aeromag). The licence is being confirmed (see
research notes). If it is the Open Government Licence: "Contains NERC materials ©NERC [year]".
