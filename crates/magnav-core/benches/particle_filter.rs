//! Particle filter throughput: particles per second for one full step (predict, weigh,
//! normalise, and resample when triggered), single-threaded vs rayon.
//!
//! Run with `cargo bench -p magnav-core`. Criterion reports throughput in elements/s, where
//! one element is one particle.

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use magnav_core::synthetic::{GaussianSumField, RandomFieldParams};
use magnav_core::{Bounds, ExecMode, FieldMap, FilterConfig, MapGrid, ParticleFilter, Vec2};

/// Same synthetic strip as the sanity tests (synthetic test values, not physical claims).
fn bench_map() -> MapGrid {
    let params = RandomFieldParams {
        bounds: Bounds {
            x_min: 0.0,
            x_max: 100_000.0,
            y_min: 40_000.0,
            y_max: 60_000.0,
        },
        n_anomalies: 200,
        amplitude_range_nt: (100.0, 500.0),
        sigma_range_m: (1_000.0, 3_000.0),
    };
    GaussianSumField::random(&params, 1).to_grid(100.0).unwrap()
}

fn filter_step(c: &mut Criterion) {
    let map = bench_map();
    let centre = Vec2::new(50_000.0, 50_000.0);
    let measurement = map.field_at(centre.x, centre.y).unwrap();

    let mut group = c.benchmark_group("filter_step");
    for n in [1_000, 10_000, 100_000] {
        group.throughput(Throughput::Elements(n as u64));
        for (label, mode) in [("serial", ExecMode::Serial), ("rayon", ExecMode::Parallel)] {
            let config = FilterConfig {
                n_particles: n,
                initial_sigma_m: 1_000.0,
                process_noise_m_per_sqrt_s: 5.0,
                measurement_sigma_nt: 2.0,
                resample_fraction: 0.5,
                seed: 1,
                mode,
            };
            let mut filter = ParticleFilter::new(config, centre);
            // Zero INS velocity keeps the cloud near the centre however many iterations
            // criterion runs, so particles never leave the map.
            group.bench_function(BenchmarkId::new(label, n), |b| {
                b.iter(|| {
                    filter
                        .step(black_box(&map), Vec2::default(), 1.0, measurement)
                        .unwrap()
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, filter_step);
criterion_main!(benches);
