//! Runs a complete simulated flight: truth → measurements → filter → per-step records.
//!
//! This is the only module that knows the true position. It hands the filter nothing but
//! the INS velocity and the magnetometer reading.

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use thiserror::Error;

use crate::filter::{FilterConfig, FilterError, ParticleFilter};
use crate::map::FieldMap;
use crate::sensor::Magnetometer;
use crate::trajectory::{ConstantBiasIns, StraightLine};
use crate::vec2::Vec2;

/// A simulated flight.
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    /// True path.
    pub truth: StraightLine,
    /// INS error model.
    pub ins: ConstantBiasIns,
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

/// Simulates `scenario` over `map` with `sensor`, running a filter configured by `config`.
///
/// The same `map` both generates measurements and is used by the filter, i.e. the map is
/// perfect. Map error is studied in Phase 2. `sensor_seed` drives the sensor noise only.
pub fn run_scenario<M, S>(
    map: &M,
    sensor: &S,
    scenario: &Scenario,
    config: &FilterConfig,
    sensor_seed: u64,
) -> Result<Vec<StepRecord>, SimError>
where
    M: FieldMap,
    S: Magnetometer,
{
    let mut sensor_rng = Xoshiro256PlusPlus::seed_from_u64(sensor_seed);
    let mut filter = ParticleFilter::new(config.clone(), scenario.ins_start);
    let ins_velocity = scenario.ins.reported_velocity(scenario.truth.velocity);

    let mut records = Vec::with_capacity(scenario.n_steps);
    for k in 1..=scenario.n_steps {
        let t = k as f64 * scenario.dt;
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
            ins: scenario
                .ins
                .dead_reckon(scenario.ins_start, &scenario.truth, t),
            estimate: summary.estimate,
            spread_m: summary.spread_m,
            ess: summary.ess,
            resampled: summary.resampled,
        });
    }
    Ok(records)
}
