//! `magnav`: command-line driver for magnav-core simulations.
//!
//! ```text
//! magnav demo --map gaussian --particles 5000 --seed 1 --out run.csv
//! magnav demo --map flat
//! magnav demo --map-file data/grid/bgs_250m.npy --start-x 200000 --start-y 750000
//! magnav demo --map-file data/grid/bgs_250m.npy --start-x 144000 --start-y 730000 \
//!     --ins hg4930 --sensor-noise-nt 2
//! ```
//!
//! Error handling: the library returns typed errors (`thiserror` enums) so callers can match
//! on them. An application like this one only needs to report them, so it uses `anyhow`.
//! Any error type converts into `anyhow::Error` via `?`, and `.context(...)` adds a message.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use magnav_core::ins::{ImuSpec, InsModel, SchulerIns};
use magnav_core::sensor::{IdealSensor, Magnetometer, WhiteNoiseSensor};
use magnav_core::synthetic::{FlatField, GaussianSumField, RandomFieldParams};
use magnav_core::trajectory::{ConstantBiasIns, StraightLine};
use magnav_core::{
    Bounds, ExecMode, FieldMap, FilterConfig, MapGrid, Scenario, StepRecord, Vec2, run_scenario,
};

#[derive(Parser)]
#[command(version, about = "Magnetic anomaly navigation simulator")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Fly a straight 60 km line east and compare the filter with the INS.
    Demo(DemoArgs),
}

#[derive(clap::Args)]
struct DemoArgs {
    /// Synthetic map type (ignored when --map-file is given).
    #[arg(long, value_enum, default_value_t = MapKind::Gaussian)]
    map: MapKind,
    /// Fly over a real gridded map (`.npy` + `.json` sidecar) instead of a synthetic one.
    #[arg(long, conflicts_with = "map", requires_all = ["start_x", "start_y"])]
    map_file: Option<PathBuf>,
    /// Start easting of the true trajectory (m, map coordinates). Default: 20000.
    #[arg(long)]
    start_x: Option<f64>,
    /// Start northing of the true trajectory (m, map coordinates). Default: 50000.
    #[arg(long)]
    start_y: Option<f64>,
    /// Number of particles.
    #[arg(long, default_value_t = 5_000)]
    particles: usize,
    /// Seed for the map, the filter and the sensor noise.
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Magnetometer white noise, standard deviation per sample (nT). 0 = ideal sensor.
    #[arg(long, default_value_t = 0.0)]
    sensor_noise_nt: f64,
    /// INS error model.
    #[arg(long, value_enum, default_value_t = InsKind::ConstantBias)]
    ins: InsKind,
    /// Filter position process noise (m/√s). Default: 5 (constant-bias), 1 (hg4930).
    #[arg(long)]
    process_noise: Option<f64>,
    /// Prior standard deviation of the INS velocity bias the filter estimates (m/s per axis).
    /// 0 together with --bias-noise 0 turns bias estimation off. Default: 2 (constant-bias),
    /// 0 (hg4930, whose velocity error is zero when GNSS is lost).
    #[arg(long)]
    bias_prior: Option<f64>,
    /// Random walk of the filter's bias model (m/s per √s). Default: 0.01 (constant-bias);
    /// for hg4930, sized from the IMU model as sigma_v(T) / sqrt(T) over the flight.
    #[arg(long)]
    bias_noise: Option<f64>,
    /// Run the filter single-threaded instead of with rayon.
    #[arg(long)]
    serial: bool,
    /// Write per-step results to this CSV file.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy, ValueEnum)]
enum InsKind {
    /// Synthetic: 2 m/s constant velocity bias, INS starts 721 m off (as in the sanity tests).
    ConstantBias,
    /// Unaided INS on a Honeywell HG4930 (CA51) IMU, Schuler-loop errors from zero at GNSS
    /// loss; the filter starts on the true position with a 10 m prior.
    Hg4930,
}

#[derive(Clone, Copy, ValueEnum)]
enum MapKind {
    /// 200 random Gaussian anomalies, sampled onto a 100 m grid.
    Gaussian,
    /// A constant field: no position information at all.
    Flat,
}

// All numbers below are synthetic demo values, not physical claims (same as the sanity tests).
const BOUNDS: Bounds = Bounds {
    x_min: 0.0,
    x_max: 100_000.0,
    y_min: 40_000.0,
    y_max: 60_000.0,
};

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Demo(args) => demo(&args),
    }
}

fn demo(args: &DemoArgs) -> Result<()> {
    if let Some(path) = &args.map_file {
        let (grid, meta) =
            MapGrid::load(path).with_context(|| format!("loading map {}", path.display()))?;
        println!(
            "map: {} ({} x {} nodes, {} m spacing)",
            meta.source.as_deref().unwrap_or("no source given"),
            grid.nx(),
            grid.ny(),
            grid.dx()
        );
        // The same map generates the measurements and drives the filter, so map error is
        // zero here: a best case. Separate truth and filter maps come in Phase 2.
        println!("note: truth and filter share one map (no map error modelled)");
        let records = simulate(&grid, args)?;
        return finish(args, &records);
    }

    let records = match args.map {
        MapKind::Gaussian => {
            let params = RandomFieldParams {
                bounds: BOUNDS,
                n_anomalies: 200,
                amplitude_range_nt: (100.0, 500.0),
                sigma_range_m: (1_000.0, 3_000.0),
            };
            let grid = GaussianSumField::random(&params, args.seed).to_grid(100.0)?;
            simulate(&grid, args)?
        }
        MapKind::Flat => {
            let flat = FlatField {
                value_nt: 50.0,
                bounds: BOUNDS,
            };
            simulate(&flat, args)?
        }
    };
    finish(args, &records)
}

/// Writes the optional CSV and prints the summary.
fn finish(args: &DemoArgs, records: &[StepRecord]) -> Result<()> {
    if let Some(path) = &args.out {
        write_csv(path, records).with_context(|| format!("writing {}", path.display()))?;
        println!("wrote {} steps to {}", records.len(), path.display());
    }
    print_summary(records);
    Ok(())
}

/// Picks the sensor model, then runs the scenario. The map type `M` and sensor type are
/// generic, so each combination gets its own compiled (and fully inlined) copy.
fn simulate<M: FieldMap + Sync>(map: &M, args: &DemoArgs) -> Result<Vec<StepRecord>> {
    if args.sensor_noise_nt > 0.0 {
        let sensor = WhiteNoiseSensor {
            sigma_nt: args.sensor_noise_nt,
        };
        run(map, &sensor, args)
    } else {
        run(map, &IdealSensor, args)
    }
}

fn run<M: FieldMap + Sync, S: Magnetometer>(
    map: &M,
    sensor: &S,
    args: &DemoArgs,
) -> Result<Vec<StepRecord>> {
    let truth = StraightLine {
        start: Vec2::new(
            args.start_x.unwrap_or(20_000.0),
            args.start_y.unwrap_or(50_000.0),
        ),
        velocity: Vec2::new(100.0, 0.0),
    };
    let (dt, n_steps) = (1.0, 600);
    // Per-INS defaults: (INS model, INS start, initial sigma, process noise, bias prior,
    // bias noise). Command-line flags override the last three.
    let (ins, ins_start, initial_sigma, q, bias_prior, bias_noise): (InsModel, _, _, _, _, _) =
        match args.ins {
            InsKind::ConstantBias => (
                ConstantBiasIns {
                    velocity_bias: Vec2::new(0.0, 2.0),
                }
                .into(),
                truth.start + Vec2::new(600.0, -400.0),
                1_000.0,
                5.0,
                2.0,
                0.01,
            ),
            InsKind::Hg4930 => {
                let ins = SchulerIns {
                    spec: ImuSpec::HG4930_CA51,
                };
                let flight_s = dt * n_steps as f64;
                let sized = ins.velocity_error_sigma(flight_s) / flight_s.sqrt();
                (InsModel::Schuler(ins), truth.start, 10.0, 1.0, 0.0, sized)
            }
        };
    let scenario = Scenario {
        truth,
        ins,
        ins_start,
        dt,
        n_steps,
    };
    let config = FilterConfig {
        n_particles: args.particles,
        initial_sigma_m: initial_sigma,
        process_noise_m_per_sqrt_s: args.process_noise.unwrap_or(q),
        initial_bias_sigma_m_per_s: args.bias_prior.unwrap_or(bias_prior),
        bias_noise_m_per_s_per_sqrt_s: args.bias_noise.unwrap_or(bias_noise),
        // The likelihood must be at least as wide as the sensor noise.
        measurement_sigma_nt: 2.0_f64.max(args.sensor_noise_nt),
        resample_fraction: 0.5,
        seed: args.seed,
        mode: if args.serial {
            ExecMode::Serial
        } else {
            ExecMode::Parallel
        },
    };
    Ok(run_scenario(map, sensor, &scenario, &config, args.seed)?)
}

fn write_csv(path: &Path, records: &[StepRecord]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "t_s,truth_x_m,truth_y_m,ins_x_m,ins_y_m,est_x_m,est_y_m,filter_error_m,ins_error_m,spread_m,nees,bias_x_m_per_s,bias_y_m_per_s,ess,resampled"
    )?;
    for r in records {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            r.t,
            r.truth.x,
            r.truth.y,
            r.ins.x,
            r.ins.y,
            r.estimate.x,
            r.estimate.y,
            r.filter_error_m(),
            r.ins_error_m(),
            r.spread_m,
            r.nees(),
            r.bias_estimate.x,
            r.bias_estimate.y,
            r.ess,
            r.resampled
        )?;
    }
    Ok(())
}

fn print_summary(records: &[StepRecord]) {
    let Some(last) = records.last() else {
        println!("no steps simulated");
        return;
    };
    let tail = &records[records.len().saturating_sub(60)..];
    let tail_mean = tail.iter().map(StepRecord::filter_error_m).sum::<f64>() / tail.len() as f64;
    println!("after {:.0} s:", last.t);
    println!("  INS-only error      {:8.1} m", last.ins_error_m());
    println!("  filter error        {:8.1} m", last.filter_error_m());
    println!("  filter, last 60 s   {:8.1} m (mean)", tail_mean);
    println!(
        "  filter spread       {:8.1} m (1-sigma per axis)",
        last.spread_m
    );
    let tail_nees = tail.iter().map(StepRecord::nees).sum::<f64>() / tail.len() as f64;
    println!("  NEES, last 60 s     {tail_nees:8.2}   (about 2 if the spread is honest)");
    println!(
        "  bias estimate       ({:+.3}, {:+.3}) m/s",
        last.bias_estimate.x, last.bias_estimate.y
    );
}
