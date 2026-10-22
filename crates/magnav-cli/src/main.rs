//! `magnav`: command-line driver for magnav-core simulations.
//!
//! ```text
//! magnav demo --map gaussian --particles 5000 --seed 1 --out run.csv
//! magnav demo --map flat
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
use magnav_core::sensor::{IdealSensor, Magnetometer, WhiteNoiseSensor};
use magnav_core::synthetic::{FlatField, GaussianSumField, RandomFieldParams};
use magnav_core::trajectory::{ConstantBiasIns, StraightLine};
use magnav_core::{
    Bounds, ExecMode, FieldMap, FilterConfig, Scenario, StepRecord, Vec2, run_scenario,
};

#[derive(Parser)]
#[command(version, about = "Magnetic anomaly navigation simulator")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Fly a straight 60 km line over a synthetic map and compare the filter with the INS.
    Demo(DemoArgs),
}

#[derive(clap::Args)]
struct DemoArgs {
    /// Synthetic map type.
    #[arg(long, value_enum, default_value_t = MapKind::Gaussian)]
    map: MapKind,
    /// Number of particles.
    #[arg(long, default_value_t = 5_000)]
    particles: usize,
    /// Seed for the map, the filter and the sensor noise.
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Magnetometer white noise, standard deviation per sample (nT). 0 = ideal sensor.
    #[arg(long, default_value_t = 0.0)]
    sensor_noise_nt: f64,
    /// Run the filter single-threaded instead of with rayon.
    #[arg(long)]
    serial: bool,
    /// Write per-step results to this CSV file.
    #[arg(long)]
    out: Option<PathBuf>,
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

    if let Some(path) = &args.out {
        write_csv(path, &records).with_context(|| format!("writing {}", path.display()))?;
        println!("wrote {} steps to {}", records.len(), path.display());
    }
    print_summary(&records);
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
        start: Vec2::new(20_000.0, 50_000.0),
        velocity: Vec2::new(100.0, 0.0),
    };
    let scenario = Scenario {
        truth,
        ins: ConstantBiasIns {
            velocity_bias: Vec2::new(0.0, 2.0),
        },
        ins_start: truth.start + Vec2::new(600.0, -400.0),
        dt: 1.0,
        n_steps: 600,
    };
    let config = FilterConfig {
        n_particles: args.particles,
        initial_sigma_m: 1_000.0,
        process_noise_m_per_sqrt_s: 5.0,
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
        "t_s,truth_x_m,truth_y_m,ins_x_m,ins_y_m,est_x_m,est_y_m,filter_error_m,ins_error_m,spread_m,ess,resampled"
    )?;
    for r in records {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{},{}",
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
}
