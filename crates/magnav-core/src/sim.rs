//! Runs a complete simulated flight: truth → measurements → filter → per-step records.
//!
//! This is the only module that knows the true position. It hands the filter nothing but
//! the INS velocity and the magnetometer reading.

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use thiserror::Error;

use crate::filter::{Cov2, FilterConfig, FilterError, ParticleFilter};
use crate::ins::InsModel;
use crate::map::FieldMap;
use crate::sensor::Magnetometer;
use crate::trajectory::StraightLine;
use crate::vec2::Vec2;

/// A simulated flight.
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    /// True path.
    pub truth: StraightLine,
    /// INS error model.
    pub ins: InsModel,
    /// Where the INS (and the filter's prior) believes the flight starts.
    pub ins_start: Vec2,
    /// Time step (s), equal to the magnetometer sample interval.
    pub dt: f64,
    /// Number of steps.
    pub n_steps: usize,
}

/// One time step of a simulated flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepRecord {
    /// Time (s).
    pub t: f64,
    /// True position.
    pub truth: Vec2,
    /// INS-only dead-reckoned position.
    pub ins: Vec2,
    /// Filter estimate.
    pub estimate: Vec2,
    /// Filter estimate of the INS velocity bias (m/s).
    pub bias_estimate: Vec2,
    /// Filter position covariance (m^2).
    pub cov: Cov2,
    /// Filter per-axis spread (m).
    pub spread_m: f64,
    /// Effective sample size after the update.
    pub ess: f64,
    /// Whether the filter resampled.
    pub resampled: bool,
}

impl StepRecord {
    /// Horizontal distance between estimate and truth (m).
    pub fn filter_error_m(&self) -> f64 {
        (self.estimate - self.truth).norm()
    }

    /// Normalised estimation error squared, `e^T P^-1 e`. Averages 2 for a filter whose
    /// reported covariance matches its actual errors; much larger means overconfident.
    pub fn nees(&self) -> f64 {
        self.cov.mahalanobis2(self.estimate - self.truth)
    }

    /// Horizontal distance between INS dead reckoning and truth (m).
    pub fn ins_error_m(&self) -> f64 {
        (self.ins - self.truth).norm()
    }
}

/// Errors from [`run_scenario`].
#[derive(Debug, Error, PartialEq)]
pub enum SimError {
    /// The true path left the map, so no measurement can be simulated.
    #[error("true position left the map at t = {t} s")]
    TruthOffMap {
        /// Time (s).
        t: f64,
    },
    /// The filter failed.
    #[error(transparent)]
    Filter(#[from] FilterError),
}

/// Mixed into the seed of the INS error stream so it is independent of the sensor noise.
const INS_STREAM: u64 = 0x9E37_79B9_7F4A_7C15;

/// Simulates `scenario` over `map` with `sensor`, running a filter configured by `config`.
///
/// The same `map` both generates measurements and is used by the filter, i.e. the map is
/// perfect. Map error is studied in Phase 2. `sensor_seed` drives the sensor noise and, through
/// a separate stream, the INS errors.
pub fn run_scenario<M, S>(
    map: &M,
    sensor: &S,
    scenario: &Scenario,
    config: &FilterConfig,
    sensor_seed: u64,
) -> Result<Vec<StepRecord>, SimError>
where
    M: FieldMap + Sync,
    S: Magnetometer,
{
    let mut sensor_rng = Xoshiro256PlusPlus::seed_from_u64(sensor_seed);
    let mut ins_rng = Xoshiro256PlusPlus::seed_from_u64(sensor_seed ^ INS_STREAM);
    let mut filter = ParticleFilter::new(config.clone(), scenario.ins_start);
    let velocity_errors = scenario
        .ins
        .velocity_errors(scenario.dt, scenario.n_steps, &mut ins_rng);
    let mut ins_position = scenario.ins_start;

    let mut records = Vec::with_capacity(scenario.n_steps);
    // `zip` pairs each step number with that step's INS velocity error.
    for (k, &velocity_error) in (1..=scenario.n_steps).zip(&velocity_errors) {
        let t = k as f64 * scenario.dt;
        let ins_velocity = scenario.truth.velocity + velocity_error;
        ins_position += ins_velocity * scenario.dt;
        let truth = scenario.truth.position_at(t);
        let true_field = map
            .field_at(truth.x, truth.y)
            .ok_or(SimError::TruthOffMap { t })?;
        let z = sensor.measure(true_field, t, &mut sensor_rng);

        // `?` converts a FilterError into SimError::Filter via the `#[from]` above.
        let summary = filter.step(map, ins_velocity, scenario.dt, z)?;
        records.push(StepRecord {
            t,
            truth,
            ins: ins_position,
            estimate: summary.estimate,
            bias_estimate: summary.bias_estimate,
            cov: summary.cov,
            spread_m: summary.spread_m,
            ess: summary.ess,
            resampled: summary.resampled,
        });
    }
    Ok(records)
}
