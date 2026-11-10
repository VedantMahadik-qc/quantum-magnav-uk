//! End-to-end sanity checks that the particle filter uses map information, and only map
//! information. Pass criteria were fixed before the tests were first run.
//!
//! All parameters here are **synthetic test values, not physical claims**.
//!
//! Integration tests live in `tests/` and can only use the crate's public API, exactly like
//! an external user would.

use magnav_core::ins::{ImuSpec, InsModel, SchulerIns};
use magnav_core::sensor::{IdealSensor, WhiteNoiseSensor};
use magnav_core::synthetic::{FlatField, GaussianSumField, RandomFieldParams};
use magnav_core::trajectory::{ConstantBiasIns, StraightLine};
use magnav_core::{
    Bounds, ExecMode, FieldMap, FilterConfig, Scenario, StepRecord, Vec2, run_scenario,
};

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
        }
        .into(),
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
        initial_bias_sigma_m_per_s: 0.0,
        bias_noise_m_per_s_per_sqrt_s: 0.0,
        measurement_sigma_nt: 2.0,
        resample_fraction: 0.5,
        seed,
        mode: ExecMode::Parallel,
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

fn run<M: FieldMap + Sync>(map: &M, seed: u64) -> Vec<StepRecord> {
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

#[test]
fn serial_and_parallel_are_bit_identical() {
    // Pass: identical step records (exact `==` on every f64) for serial and for rayon
    // pools of 1, 3 and 16 threads.
    let seed = 7;
    let grid = gaussian_map(seed).to_grid(100.0).unwrap();
    let with_mode = |mode| {
        let cfg = FilterConfig {
            mode,
            ..config(seed)
        };
        run_scenario(&grid, &IdealSensor, &scenario(), &cfg, seed).unwrap()
    };
    let serial = with_mode(ExecMode::Serial);
    for threads in [1, 3, 16] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let parallel = pool.install(|| with_mode(ExecMode::Parallel));
        assert_eq!(serial, parallel, "differs with {threads} threads");
    }
}

// ---- Velocity-bias estimation. Criteria fixed before the first run. ----
//
// History: the first implementation sampled the bias as extra particle coordinates. It passed
// the position and NEES checks but failed the bias check (errors 0.35-1.0 m/s while the
// filter claimed ~0.2 m/s; with less position jitter some seeds diverged to ~300 m). That
// is sample impoverishment. The filter now marginalises the bias (see `filter.rs`) and must
// pass the same criteria.

/// The true INS velocity bias in `scenario()`.
const TRUE_BIAS: Vec2 = Vec2::new(0.0, 2.0);

/// `config` with bias estimation on: a 2 m/s per-axis prior (so the true bias sits at 1 sigma)
/// and a 0.01 m/s/√s random walk in the filter's bias model. Synthetic test values.
fn bias_config(seed: u64) -> FilterConfig {
    FilterConfig {
        initial_bias_sigma_m_per_s: 2.0,
        bias_noise_m_per_s_per_sqrt_s: 0.01,
        ..config(seed)
    }
}

fn run_with<M: FieldMap + Sync>(map: &M, cfg: &FilterConfig, seed: u64) -> Vec<StepRecord> {
    run_scenario(map, &IdealSensor, &scenario(), cfg, seed).unwrap()
}

/// Mean NEES over the last 300 steps of a run.
fn tail_nees(records: &[StepRecord]) -> f64 {
    let tail = &records[records.len() - 300..];
    tail.iter().map(StepRecord::nees).sum::<f64>() / tail.len() as f64
}

#[test]
fn estimating_the_bias_beats_the_position_only_filter() {
    // Pass, every seed: final error < 0.5 x the 2-D filter's, and the final bias estimate
    // within 0.2 m/s of the truth.
    for seed in SEEDS {
        let grid = gaussian_map(seed).to_grid(100.0).unwrap();
        let without = run_with(&grid, &config(seed), seed);
        let with = run_with(&grid, &bias_config(seed), seed);
        let (e2, e4) = (
            without.last().unwrap().filter_error_m(),
            with.last().unwrap().filter_error_m(),
        );
        let bias_err = (with.last().unwrap().bias_estimate - TRUE_BIAS).norm();
        println!("seed {seed}: 2-D {e2:.1} m, with bias {e4:.1} m, bias error {bias_err:.3} m/s");
        assert!(e4 < 0.5 * e2, "seed {seed}: {e4:.1} m vs 2-D {e2:.1} m");
        assert!(bias_err < 0.2, "seed {seed}: bias error {bias_err:.3} m/s");
    }
}

#[test]
fn bias_filter_is_consistent_on_an_informative_map() {
    // Pass: mean NEES over the last 300 s, averaged over seeds, <= 6. A consistent 2-D
    // filter averages 2; the bound is loose because successive NEES values are correlated.
    let mut with = 0.0;
    let mut without = 0.0;
    for seed in SEEDS {
        let grid = gaussian_map(seed).to_grid(100.0).unwrap();
        with += tail_nees(&run_with(&grid, &bias_config(seed), seed)) / SEEDS.len() as f64;
        without += tail_nees(&run_with(&grid, &config(seed), seed)) / SEEDS.len() as f64;
    }
    println!("mean NEES, last 300 s: with bias {with:.2}, position-only {without:.2}");
    assert!(with <= 6.0, "mean NEES {with:.2}");
}

#[test]
fn bias_filter_does_not_converge_on_a_flat_map() {
    // With no map information the bias cannot be learned either. Pass, every seed: final
    // error >= 0.8 x INS-only error.
    let flat = FlatField {
        value_nt: 50.0,
        bounds: BOUNDS,
    };
    for seed in SEEDS {
        let last = *run_with(&flat, &bias_config(seed), seed).last().unwrap();
        println!(
            "flat seed {seed}: filter {:.1} m, INS {:.1} m",
            last.filter_error_m(),
            last.ins_error_m()
        );
        assert!(last.filter_error_m() >= 0.8 * last.ins_error_m());
    }
}

// ---- Realistic INS: HG4930 (CA51) free-inertial errors. ----
//
// GNSS is lost at t = 0 with the INS and filter exactly on the true position; the filter's
// prior is 10 m. The filter's bias random walk is sized from the INS model so that over the
// flight it spans the expected 1-sigma velocity error: q_b = sigma_v(T) / sqrt(T).
//
// History: the criteria below were first run with an *ideal* sensor and failed. Mean NEES was
// 0.17-0.19 for every q in the sweep, so q was not the cause. A diagnostic then showed that
// the filter's 2 nT likelihood is the problem when the readings are in fact exact: with 2 nT
// white sensor noise (matching the likelihood) NEES rose to 1.2-1.3. The test now uses
// matched noise. That change was made AFTER seeing the diagnostic, so this test is a
// regression check, not a pre-registered prediction.

const HG4930: SchulerIns = SchulerIns {
    spec: ImuSpec::HG4930_CA51,
};

fn hg4930_scenario() -> Scenario {
    let truth = scenario().truth;
    Scenario {
        ins: InsModel::Schuler(HG4930),
        ins_start: truth.start,
        ..scenario()
    }
}

fn hg4930_config(seed: u64, process_noise: f64) -> FilterConfig {
    let flight_s = 600.0;
    FilterConfig {
        initial_sigma_m: 10.0,
        process_noise_m_per_sqrt_s: process_noise,
        initial_bias_sigma_m_per_s: 0.0,
        bias_noise_m_per_s_per_sqrt_s: HG4930.velocity_error_sigma(flight_s) / flight_s.sqrt(),
        ..config(seed)
    }
}

/// One point of the process-noise sweep.
struct SweepPoint {
    q: f64,
    nees: f64,
    /// (filter, INS-only) final errors per seed (m).
    finals: Vec<(f64, f64)>,
}

#[test]
fn process_noise_sized_for_a_real_imu_gives_a_consistent_filter() {
    // Sweep q over {0.5, 1, 2, 5} m/√s. Pick the q whose mean NEES (last 300 s, averaged
    // over seeds) is closest to 2 in log ratio. Pass: that NEES is in [0.5, 8], and with it
    // every seed ends with error < 0.2 x the INS-only error.
    let mut best: Option<SweepPoint> = None;
    for q in [0.5, 1.0, 2.0, 5.0] {
        let mut nees = 0.0;
        let mut finals = Vec::new();
        for seed in SEEDS {
            let grid = gaussian_map(seed).to_grid(100.0).unwrap();
            let records = run_scenario(
                &grid,
                &WhiteNoiseSensor { sigma_nt: 2.0 },
                &hg4930_scenario(),
                &hg4930_config(seed, q),
                seed,
            )
            .unwrap();
            nees += tail_nees(&records) / SEEDS.len() as f64;
            let last = records.last().unwrap();
            finals.push((last.filter_error_m(), last.ins_error_m()));
        }
        println!("q = {q} m/√s: mean NEES {nees:.2}, (filter, INS) final errors {finals:.1?}");
        let better = best
            .as_ref()
            .is_none_or(|b| (nees / 2.0).ln().abs() < (b.nees / 2.0).ln().abs());
        if better {
            best = Some(SweepPoint { q, nees, finals });
        }
    }
    let SweepPoint { q, nees, finals } = best.unwrap();
    println!("chosen q = {q} m/√s (NEES {nees:.2})");
    assert!(
        (0.5..=8.0).contains(&nees),
        "best NEES {nees:.2} at q = {q}"
    );
    for (seed, (err, ins_err)) in SEEDS.iter().zip(finals) {
        assert!(
            err < 0.2 * ins_err,
            "seed {seed}: {err:.1} m vs INS {ins_err:.1} m"
        );
    }
}
