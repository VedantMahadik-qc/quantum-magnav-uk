//! Magnetometer measurement models.
//!
//! Phase 1 only needs an ideal sensor and white Gaussian noise. The fluxgate,
//! optically pumped and NV models (drift, heading error, calibration residuals) come in
//! Phase 2, parameterised from the sourced table in `docs/research-notes.md`.

use rand::{Rng, RngExt};
use rand_distr::StandardNormal;

/// A scalar (total-field) magnetometer.
///
/// `measure` is generic over the random number generator `R`. `?Sized` also allows
/// `&mut dyn Rng` trait objects, not just concrete generator types.
pub trait Magnetometer {
    /// One reading (nT) when the true anomaly field is `true_field_nt` at time `t` (s).
    fn measure<R: Rng + ?Sized>(&self, true_field_nt: f64, t: f64, rng: &mut R) -> f64;
}

/// A perfect sensor: returns the true field.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct IdealSensor;

impl Magnetometer for IdealSensor {
    fn measure<R: Rng + ?Sized>(&self, true_field_nt: f64, _t: f64, _rng: &mut R) -> f64 {
        true_field_nt
    }
}

/// Adds independent zero-mean Gaussian noise with standard deviation `sigma_nt` per sample.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhiteNoiseSensor {
    /// Per-sample noise standard deviation (nT).
    pub sigma_nt: f64,
}

impl WhiteNoiseSensor {
    /// Per-sample noise from a white one-sided amplitude spectral density (nT/√Hz) sampled
    /// at `sample_rate_hz`: `sigma = ASD * sqrt(f_s / 2)`. That is the noise power integrated
    /// over the band from 0 to the Nyquist frequency `f_s / 2`. For flat (white) noise, the
    /// variance is the one-sided PSD (ASD²) integrated over that band: `ASD² * f_s / 2`.
    pub fn from_asd(asd_nt_per_sqrt_hz: f64, sample_rate_hz: f64) -> Self {
        Self {
            sigma_nt: asd_nt_per_sqrt_hz * (sample_rate_hz / 2.0).sqrt(),
        }
    }
}

impl Magnetometer for WhiteNoiseSensor {
    fn measure<R: Rng + ?Sized>(&self, true_field_nt: f64, _t: f64, rng: &mut R) -> f64 {
        let z: f64 = rng.sample(StandardNormal);
        true_field_nt + self.sigma_nt * z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn ideal_sensor_is_exact() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        assert_eq!(IdealSensor.measure(123.4, 0.0, &mut rng), 123.4);
    }

    #[test]
    fn white_noise_has_requested_mean_and_std() {
        // Pass criteria: with N = 100 000 samples, the mean is within 4 standard errors
        // (4 sigma / sqrt(N)) and the std within 1% (its standard error is ~0.22%).
        let sensor = WhiteNoiseSensor { sigma_nt: 5.0 };
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let n = 100_000;
        let samples: Vec<f64> = (0..n)
            .map(|_| sensor.measure(10.0, 0.0, &mut rng))
            .collect();
        let mean = samples.iter().sum::<f64>() / n as f64;
        let var = samples.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        assert!(
            (mean - 10.0).abs() < 4.0 * 5.0 / (n as f64).sqrt(),
            "mean {mean}"
        );
        assert!((var.sqrt() / 5.0 - 1.0).abs() < 0.01, "std {}", var.sqrt());
    }

    #[test]
    fn asd_conversion_matches_nyquist_bandwidth() {
        // 1 nT/√Hz sampled at 2 Hz: band = 1 Hz, so sigma = 1 nT.
        assert_eq!(WhiteNoiseSensor::from_asd(1.0, 2.0).sigma_nt, 1.0);
        // 4 pT/√Hz (G-823A spec, research-notes) at 10 Hz -> 4e-3 * sqrt(5) nT.
        let s = WhiteNoiseSensor::from_asd(4e-3, 10.0).sigma_nt;
        assert!((s - 4e-3 * 5f64.sqrt()).abs() < 1e-15);
    }
}
