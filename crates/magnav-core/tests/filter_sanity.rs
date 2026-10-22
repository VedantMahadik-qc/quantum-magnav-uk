//! End-to-end sanity checks that the particle filter uses map information, and only map
//! information. Pass criteria were fixed before the tests were first run.
//!
//! All parameters here are **synthetic test values, not physical claims**.
//!
//! Integration tests live in `tests/` and can only use the crate's public API, exactly like
//! an external user would.

use magnav_core::sensor::IdealSensor;
use magnav_core::synthetic::{FlatField, GaussianSumField, RandomFieldParams};
use magnav_core::trajectory::{ConstantBiasIns, StraightLine};
use magnav_core::{Bounds, FieldMap, FilterConfig, Scenario, StepRecord, Vec2, run_scenario};

const SEEDS: [u64; 5] = [1, 2, 3, 4, 5];
const N_PARTICLES: usize = 5_000;
const INITIAL_SIGMA_M: f64 = 1_000.0;
const PROCESS_NOISE: f64 = 5.0; // m/√s

/// A 100 km x 20 km strip around the flight line. The filter's particles stay within a few
/// km of the line, so they never reach the edges.
const BOUNDS: Bounds = Bounds {
    x_min: 0.0,
    x_max: 100_000.0,
    y_min: 40_000.0,
    y_max: 60_000.0,
};

/// 600 s east at 100 m/s through the middle of the map. The INS starts 721 m off and has a
/// 2 m/s northward velocity bias, so its final error is |(600, -400 + 1200)| = 1000 m.
fn scenario() -> Scenario {
    let truth = StraightLine {
        start: Vec2::new(20_000.0, 50_000.0),
        velocity: Vec2::new(100.0, 0.0),
    };
    Scenario {
        truth,
        ins: ConstantBiasIns {
            velocity_bias: Vec2::new(0.0, 2.0),
        },
        ins_start: truth.start + Vec2::new(600.0, -400.0),
        dt: 1.0,
        n_steps: 600,
    }
}

fn config(seed: u64) -> FilterConfig {
    FilterConfig {
        n_particles: N_PARTICLES,
        initial_sigma_m: INITIAL_SIGMA_M,
        process_noise_m_per_sqrt_s: PROCESS_NOISE,
        measurement_sigma_nt: 2.0,
        resample_fraction: 0.5,
        seed,
    }
}

/// A strongly varying map: 200 anomalies over 2000 km^2 (mean spacing ~3 km), so the
/// track is never far from one. `map_varies_strongly_along_the_track` checks this.
fn gaussian_map(seed: u64) -> GaussianSumField {
    let params = RandomFieldParams {
        bounds: BOUNDS,
        n_anomalies: 200,
        amplitude_range_nt: (100.0, 500.0),
        sigma_range_m: (1_000.0, 3_000.0),
    };
    GaussianSumField::random(&params, seed)
}

fn run<M: FieldMap>(map: &M, seed: u64) -> Vec<StepRecord> {
    run_scenario(map, &IdealSensor, &scenario(), &config(seed), seed).unwrap()
}

#[test]
fn ins_alone_drifts_by_the_expected_amount() {
    // Checks the test set-up itself before trusting the filter comparisons below.
    let records = run(&gaussian_map(1), 1);
    let last = records.last().unwrap();
    assert!((last.ins_error_m() - 1_000.0).abs() < 1e-6);
}

#[test]
fn map_varies_strongly_along_the_track() {
    // Precondition for the convergence tests: the map must actually carry information along
    // the flight. |grad M| >= 0.01 nT/m at >= 95% of the track points.
    let path = scenario().truth;
    for seed in SEEDS {
        let map = gaussian_map(seed);
        let h = 10.0; // central-difference step (m)
        let informative = (1..=600)
            .filter(|&k| {
                let p = path.position_at(k as f64);
                let gx = (map.value_at(p.x + h, p.y) - map.value_at(p.x - h, p.y)) / (2.0 * h);
                let gy = (map.value_at(p.x, p.y + h) - map.value_at(p.x, p.y - h)) / (2.0 * h);
                gx.hypot(gy) >= 0.01
            })
            .count();
        assert!(
            informative >= 570,
            "seed {seed}: only {informative}/600 informative points"
        );
    }
}

/// Pass: final error < 100 m AND < 0.2 x INS-only final error.
fn assert_converged(label: &str, last: &StepRecord) {
    let (err, ins_err) = (last.filter_error_m(), last.ins_error_m());
    println!("{label}: filter {err:.1} m, INS {ins_err:.1} m");
    assert!(err < 100.0, "{label}: final error {err:.1} m");
    assert!(
        err < 0.2 * ins_err,
        "{label}: {err:.1} m vs INS {ins_err:.1} m"
    );
}

#[test]
fn converges_on_an_informative_map() {
    // Every seed, through MapGrid + bilinear interpolation (100 m spacing). The same grid
    // generates the measurements, so the map is perfect here; map error is Phase 2.
    for seed in SEEDS {
        let grid = gaussian_map(seed).to_grid(100.0).unwrap();
        assert_converged(
            &format!("grid seed {seed}"),
            run(&grid, seed).last().unwrap(),
        );
    }
}

#[test]
fn converges_on_the_exact_analytic_map() {
    // Same criteria on the analytic field (no interpolation at all), for one seed: it costs
    // 200 exp() calls per particle per step, so it is much slower than the grid.
    let seed = 1;
    assert_converged(
        "analytic seed 1",
        run(&gaussian_map(seed), seed).last().unwrap(),
    );
}

#[test]
fn does_not_converge_on_a_flat_map() {
    // A constant field carries no position information, so the filter must behave like pure
    // dead reckoning: it follows the INS, and its spread grows as sqrt(s0^2 + q^2 t).
    let flat = FlatField {
        value_nt: 50.0,
        bounds: BOUNDS,
    };
    let n = N_PARTICLES as f64;
    for seed in SEEDS {
        let records = run(&flat, seed);
        for r in &records {
            // (a) All likelihoods are equal, so weights stay uniform: ESS = N, no resampling.
            assert!(!r.resampled, "seed {seed}: resampled at t = {}", r.t);
            assert!((r.ess - n).abs() < 1e-9 * n, "seed {seed}: ESS {}", r.ess);

            // (d) Spread matches the analytic dead-reckoning spread to within 5%.
            let expected = (INITIAL_SIGMA_M.powi(2) + PROCESS_NOISE.powi(2) * r.t).sqrt();
            assert!(
                (r.spread_m / expected - 1.0).abs() < 0.05,
                "seed {seed}, t = {}: spread {:.1} m, expected {expected:.1} m",
                r.t,
                r.spread_m
            );

            // (b) The estimate tracks INS dead reckoning within 4 standard errors per axis.
            let tolerance = 4.0 * expected / n.sqrt();
            let offset = r.estimate - r.ins;
            assert!(
                offset.x.abs() < tolerance && offset.y.abs() < tolerance,
                "seed {seed}, t = {}: estimate is {offset:?} from INS (tolerance {tolerance:.1} m)",
                r.t
            );
        }
        // (c) The final error stays at least 80% of the INS-only error: no convergence.
        let last = records.last().unwrap();
        println!(
            "flat seed {seed}: filter {:.1} m, INS {:.1} m",
            last.filter_error_m(),
            last.ins_error_m()
        );
        assert!(last.filter_error_m() >= 0.8 * last.ins_error_m());
    }
}
