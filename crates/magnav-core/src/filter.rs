//! A bootstrap particle filter for 2-D map-matching navigation.
//!
//! Each particle is a hypothesis about the vehicle's horizontal position. Every step:
//!
//! 1. **Predict:** move each particle by the INS-reported velocity times `dt`, plus random
//!    jitter (process noise) so the cloud can absorb INS errors.
//! 2. **Update:** multiply each particle's weight by the likelihood of the magnetometer
//!    reading given the map value at that particle: `exp(-(z - M(x))^2 / (2 sigma_R^2))`.
//!    Computed in log space, subtracting the maximum before `exp`, to avoid underflow.
//! 3. **Resample** (systematic; see Douc & Cappé, "Comparison of resampling schemes for
//!    particle filtering", ISPA 2005) when the effective sample size drops below a
//!    threshold. This replaces low-weight particles with copies of high-weight ones.
//!
//! The filter never sees the true position. Its only inputs are the map, the INS velocity
//! and the measurement.
//!
//! ## Reproducibility
//! Particle `i` at step `k` draws its random numbers from its own generator seeded by
//! `(seed, k, i)`. Results therefore do not depend on the order particles are processed in,
//! which is what lets the parallel mode give bit-identical answers to the serial one.

use rand::{RngExt, SeedableRng};
use rand_distr::StandardNormal;
use rand_xoshiro::Xoshiro256PlusPlus;
use thiserror::Error;

use crate::map::FieldMap;
use crate::vec2::Vec2;

/// Tuning parameters for [`ParticleFilter`].
#[derive(Debug, Clone, PartialEq)]
pub struct FilterConfig {
    /// Number of particles.
    pub n_particles: usize,
    /// Standard deviation (m, per axis) of the initial position uncertainty.
    pub initial_sigma_m: f64,
    /// Process noise: per-axis position random walk, in m/√s. A step of `dt` seconds adds
    /// jitter with standard deviation `process_noise * sqrt(dt)`.
    pub process_noise_m_per_sqrt_s: f64,
    /// Measurement likelihood width (nT). It should cover sensor noise plus map error.
    pub measurement_sigma_nt: f64,
    /// Resample when ESS < `resample_fraction * n_particles`.
    pub resample_fraction: f64,
    /// Seed for all of the filter's randomness.
    pub seed: u64,
}

/// Errors the filter can hit while running.
#[derive(Debug, Error, PartialEq)]
pub enum FilterError {
    /// Every particle is off the map, so there is nothing left to weight. The filter state
    /// is not usable after this.
    #[error("all particles left the map at step {step}")]
    AllParticlesOffMap {
        /// Step number (1-based) at which it happened.
        step: u64,
    },
}

/// What happened in one filter step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepSummary {
    /// Weighted mean position (m).
    pub estimate: Vec2,
    /// Per-axis RMS spread of the weighted particle cloud (m): `sqrt(trace(cov) / 2)`.
    pub spread_m: f64,
    /// Effective sample size `1 / sum(w^2)` after the update (before any resampling).
    pub ess: f64,
    /// Whether resampling happened this step.
    pub resampled: bool,
}

/// The particle filter state.
#[derive(Debug, Clone)]
pub struct ParticleFilter {
    config: FilterConfig,
    positions: Vec<Vec2>,
    /// Normalised weights, summing to 1.
    weights: Vec<f64>,
    /// Number of completed steps.
    step: u64,
    /// Scratch buffers reused every step, to avoid allocating in the hot loop.
    log_weights: Vec<f64>,
    resample_indices: Vec<usize>,
}

/// Index used for the resampling random draw, so it never collides with a particle index.
const RESAMPLE_STREAM: u64 = u64::MAX;

impl ParticleFilter {
    /// Creates a filter with particles drawn from `N(initial_guess, initial_sigma^2 I)`.
    pub fn new(config: FilterConfig, initial_guess: Vec2) -> Self {
        let n = config.n_particles;
        assert!(n > 0, "need at least one particle");
        let positions = (0..n)
            .map(|i| {
                let mut rng = particle_rng(config.seed, 0, i as u64);
                initial_guess + gaussian_vec2(&mut rng) * config.initial_sigma_m
            })
            .collect();
        Self {
            weights: vec![1.0 / n as f64; n],
            log_weights: vec![0.0; n],
            resample_indices: Vec::with_capacity(n),
            positions,
            step: 0,
            config,
        }
    }

    /// Particle positions (m).
    pub fn positions(&self) -> &[Vec2] {
        &self.positions
    }

    /// Normalised particle weights.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// The filter's configuration.
    pub fn config(&self) -> &FilterConfig {
        &self.config
    }

    /// Weighted mean position.
    pub fn estimate(&self) -> Vec2 {
        self.positions
            .iter()
            .zip(&self.weights)
            .fold(Vec2::default(), |acc, (&p, &w)| acc + p * w)
    }

    /// Advances the filter by one time step of `dt` seconds.
    ///
    /// `map` is borrowed (`&M`): the filter reads it but never owns or changes it.
    /// `&mut self` means this method has exclusive access to the filter while it runs.
    pub fn step<M: FieldMap>(
        &mut self,
        map: &M,
        ins_velocity: Vec2,
        dt: f64,
        measurement_nt: f64,
    ) -> Result<StepSummary, FilterError> {
        self.step += 1;
        let step = self.step;
        let seed = self.config.seed;
        let jitter_sigma = self.config.process_noise_m_per_sqrt_s * dt.sqrt();
        let inv_two_r2 = 1.0 / (2.0 * self.config.measurement_sigma_nt.powi(2));
        let displacement = ins_velocity * dt;

        // Predict + update, one particle at a time. Each particle only touches its own
        // position, weight and random stream.
        for (i, ((pos, &w), log_w)) in self
            .positions
            .iter_mut()
            .zip(&self.weights)
            .zip(&mut self.log_weights)
            .enumerate()
        {
            let mut rng = particle_rng(seed, step, i as u64);
            *pos += displacement + gaussian_vec2(&mut rng) * jitter_sigma;
            *log_w = match map.field_at(pos.x, pos.y) {
                Some(predicted) => w.ln() - (measurement_nt - predicted).powi(2) * inv_two_r2,
                None => f64::NEG_INFINITY, // off the map: impossible
            };
        }

        // Normalise: w_i = exp(lw_i - max) / sum. Subtracting the max keeps the largest
        // term at exp(0) = 1, so we never divide 0 by 0 even for huge mismatches.
        let max_log_w = self
            .log_weights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        if max_log_w == f64::NEG_INFINITY {
            return Err(FilterError::AllParticlesOffMap { step });
        }
        for (w, &lw) in self.weights.iter_mut().zip(&self.log_weights) {
            *w = (lw - max_log_w).exp();
        }
        let total: f64 = self.weights.iter().sum();
        self.weights.iter_mut().for_each(|w| *w /= total);

        let ess = 1.0 / self.weights.iter().map(|w| w * w).sum::<f64>();
        let estimate = self.estimate();
        let spread_m = self.spread(estimate);

        let resampled = ess < self.config.resample_fraction * self.config.n_particles as f64;
        if resampled {
            let u0: f64 = particle_rng(seed, step, RESAMPLE_STREAM).random();
            systematic_resample(&self.weights, u0, &mut self.resample_indices);
            // Gather the chosen particles. `resample_indices` is a separate buffer, so we can
            // read from `positions` while building the new list.
            let new_positions: Vec<Vec2> = self
                .resample_indices
                .iter()
                .map(|&j| self.positions[j])
                .collect();
            self.positions = new_positions;
            self.weights.fill(1.0 / self.config.n_particles as f64);
        }

        Ok(StepSummary {
            estimate,
            spread_m,
            ess,
            resampled,
        })
    }

    /// Per-axis RMS spread of the weighted cloud around `mean`.
    fn spread(&self, mean: Vec2) -> f64 {
        let trace: f64 = self
            .positions
            .iter()
            .zip(&self.weights)
            .map(|(&p, &w)| {
                let d = p - mean;
                w * (d.x * d.x + d.y * d.y)
            })
            .sum();
        (trace / 2.0).sqrt()
    }
}

/// Two independent standard-normal draws as a vector.
fn gaussian_vec2(rng: &mut Xoshiro256PlusPlus) -> Vec2 {
    Vec2::new(rng.sample(StandardNormal), rng.sample(StandardNormal))
}

/// SplitMix64 finaliser (Steele, Lea & Flood, "Fast splittable pseudorandom number
/// generators", OOPSLA 2014). It scrambles a 64-bit integer so nearby inputs give unrelated
/// outputs.
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Independent random stream for particle `index` at step `step`.
fn particle_rng(seed: u64, step: u64, index: u64) -> Xoshiro256PlusPlus {
    let key = splitmix64(seed ^ splitmix64(step ^ splitmix64(index)));
    Xoshiro256PlusPlus::seed_from_u64(key)
}

/// Systematic resampling: one uniform offset `u0` in [0, 1) and N evenly spaced pointers
/// `(u0 + i) / N` walked along the cumulative weights. Low variance and O(N).
/// Particles with zero weight are never selected.
fn systematic_resample(weights: &[f64], u0: f64, out: &mut Vec<usize>) {
    let n = weights.len();
    out.clear();
    // Use the actual total rather than assuming exactly 1.0, so rounding can't push a
    // pointer past the end.
    let total: f64 = weights.iter().sum();
    let mut cumulative = weights[0];
    let mut j = 0;
    for i in 0..n {
        let u = (u0 + i as f64) / n as f64 * total;
        while cumulative <= u {
            j += 1;
            cumulative += weights[j];
        }
        out.push(j);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::map::Bounds;
    use crate::synthetic::FlatField;

    fn config(n: usize) -> FilterConfig {
        FilterConfig {
            n_particles: n,
            initial_sigma_m: 100.0,
            process_noise_m_per_sqrt_s: 1.0,
            measurement_sigma_nt: 2.0,
            resample_fraction: 0.5,
            seed: 1,
        }
    }

    #[test]
    fn systematic_resampling_picks_only_weighted_particles() {
        // Weights [0.5, 0.5, 0, 0] must give each of particles 0 and 1 exactly N/2 copies,
        // whatever the random offset.
        let weights = [0.5, 0.5, 0.0, 0.0];
        let mut out = Vec::new();
        for u0 in [0.0, 0.1, 0.5, 0.999_999] {
            systematic_resample(&weights, u0, &mut out);
            assert_eq!(out, vec![0, 0, 1, 1], "u0 = {u0}");
        }
    }

    #[test]
    fn systematic_resampling_matches_expected_counts() {
        // Each particle gets floor(N w) or ceil(N w) copies, by construction.
        let weights = [0.1, 0.0, 0.25, 0.65];
        let mut out = Vec::new();
        systematic_resample(&weights, 0.37, &mut out);
        assert_eq!(out.len(), 4);
        for (k, &w) in weights.iter().enumerate() {
            let count = out.iter().filter(|&&j| j == k).count() as f64;
            assert!(
                (count - 4.0 * w).abs() < 1.0,
                "particle {k}: {count} copies"
            );
        }
    }

    #[test]
    fn weights_stay_normalised() {
        let bounds = Bounds {
            x_min: -1e5,
            x_max: 1e5,
            y_min: -1e5,
            y_max: 1e5,
        };
        // A map with a gradient so the weights are genuinely non-uniform.
        let map =
            crate::MapGrid::from_fn(-1e5, -1e5, 1e3, 1e3, 201, 201, |x, y| 0.01 * x + 0.02 * y)
                .unwrap();
        assert_eq!(map.bounds(), bounds);
        let mut pf = ParticleFilter::new(config(1000), Vec2::default());
        for _ in 0..5 {
            pf.step(&map, Vec2::new(10.0, 0.0), 1.0, 0.0).unwrap();
            let sum: f64 = pf.weights().iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "weights sum to {sum}");
        }
    }

    #[test]
    fn errors_when_every_particle_is_off_the_map() {
        let tiny = FlatField {
            value_nt: 0.0,
            bounds: Bounds {
                x_min: 1e6,
                x_max: 1e6 + 1.0,
                y_min: 0.0,
                y_max: 1.0,
            },
        };
        let mut pf = ParticleFilter::new(config(100), Vec2::default());
        assert_eq!(
            pf.step(&tiny, Vec2::default(), 1.0, 0.0),
            Err(FilterError::AllParticlesOffMap { step: 1 })
        );
    }

    #[test]
    fn particle_streams_are_distinct() {
        let a: u64 = particle_rng(1, 1, 0).random();
        let b: u64 = particle_rng(1, 1, 1).random();
        let c: u64 = particle_rng(1, 2, 0).random();
        let d: u64 = particle_rng(2, 1, 0).random();
        assert!(a != b && a != c && a != d && b != c);
        assert_eq!(a, particle_rng(1, 1, 0).random::<u64>());
    }
}
