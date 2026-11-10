//! Free-inertial (unaided) INS error growth: the single-axis Schuler loop.
//!
//! After GNSS is lost, an INS integrates its accelerometers and gyros alone. Each horizontal
//! axis then has a velocity error `dv` and a tilt error `dth` obeying (a standard textbook
//! model, e.g. Groves, *Principles of GNSS, Inertial, and Multisensor Integrated Navigation
//! Systems*, 2nd ed., 2013):
//!
//! ```text
//! d(dv)/dt  = -g * dth + b_a + w_a      accelerometer bias b_a, white noise w_a (VRW)
//! d(dth)/dt =  dv / R  + eps + w_g      gyro bias eps, white noise w_g (ARW)
//! ```
//!
//! A tilt `dth` makes the INS see a fraction `g * dth` of gravity as horizontal acceleration.
//! The `dv / R` term (R = Earth radius) is what bounds the error into an oscillation with the
//! Schuler period `2 pi sqrt(R / g)` (about 84 min).
//!
//! Simplifications, all acceptable for flights much shorter than a Schuler period: the two
//! axes are independent (no Earth-rate or Coriolis cross-coupling), there is no vertical
//! channel, and the errors are zero at the moment GNSS is lost. The biases are drawn once per
//! flight at the IMU's **in-run bias stability**, i.e. what is left after a GNSS-aided INS has
//! calibrated out the larger turn-on bias.

use rand::{Rng, RngExt};
use rand_distr::StandardNormal;

use crate::trajectory::ConstantBiasIns;
use crate::vec2::Vec2;

/// Standard gravity (m/s^2).
pub const GRAVITY_M_S2: f64 = 9.806_65;
/// Mean Earth radius (m).
pub const EARTH_RADIUS_M: f64 = 6_371_000.0;

/// Schuler angular frequency `sqrt(g / R)` (rad/s).
pub fn schuler_omega() -> f64 {
    (GRAVITY_M_S2 / EARTH_RADIUS_M).sqrt()
}

/// IMU error parameters, in the units datasheets use.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImuSpec {
    /// Gyro in-run bias stability (deg/h, 1 sigma).
    pub gyro_bias_deg_per_h: f64,
    /// Gyro angle random walk (deg/√h).
    pub arw_deg_per_sqrt_h: f64,
    /// Accelerometer in-run bias stability (mg, 1 sigma; 1 mg = 1e-3 * 9.80665 m/s^2).
    pub accel_bias_mg: f64,
    /// Accelerometer velocity random walk (m/s/√h).
    pub vrw_m_per_s_per_sqrt_h: f64,
}

impl ImuSpec {
    /// Honeywell HG4930 tactical-grade MEMS IMU, grade CA51 (the third column of the tables
    /// in *HG4930 Inertial Measurement Unit: Performance and Environmental Information*,
    /// N61-1775-000-000, 06/17; Tables 1, 2, 4, 5). Bias values are typical 1 sigma; ARW and
    /// VRW are the stated maxima.
    pub const HG4930_CA51: ImuSpec = ImuSpec {
        gyro_bias_deg_per_h: 1.0,
        arw_deg_per_sqrt_h: 0.09,
        accel_bias_mg: 0.3,
        vrw_m_per_s_per_sqrt_h: 0.09,
    };

    fn gyro_bias_rad_s(&self) -> f64 {
        self.gyro_bias_deg_per_h.to_radians() / 3600.0
    }
    fn arw_rad_per_sqrt_s(&self) -> f64 {
        self.arw_deg_per_sqrt_h.to_radians() / 60.0
    }
    fn accel_bias_m_s2(&self) -> f64 {
        self.accel_bias_mg * 1e-3 * GRAVITY_M_S2
    }
    fn vrw_m_per_s_per_sqrt_s(&self) -> f64 {
        self.vrw_m_per_s_per_sqrt_h / 60.0
    }
}

/// The biases of one simulated flight, per axis (x = east, y = north).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InsBiases {
    /// Accelerometer bias (m/s^2).
    pub accel_m_s2: Vec2,
    /// Gyro bias (rad/s).
    pub gyro_rad_s: Vec2,
}

/// An unaided INS built on a given IMU.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SchulerIns {
    /// The IMU's error parameters.
    pub spec: ImuSpec,
}

impl SchulerIns {
    /// Draws the biases of one flight: independent zero-mean Gaussians per axis with the
    /// spec's 1 sigma values.
    pub fn draw_biases<R: Rng + ?Sized>(&self, rng: &mut R) -> InsBiases {
        let mut g = || -> f64 { rng.sample(StandardNormal) };
        let (sa, sg) = (self.spec.accel_bias_m_s2(), self.spec.gyro_bias_rad_s());
        InsBiases {
            accel_m_s2: Vec2::new(g() * sa, g() * sa),
            gyro_rad_s: Vec2::new(g() * sg, g() * sg),
        }
    }

    /// Velocity errors at `t = dt, 2 dt, ..., n_steps dt`, for given biases, starting from zero
    /// error. White noise (VRW, ARW) is drawn from `rng`; pass `with_noise = false` for the
    /// deterministic bias-only response.
    ///
    /// Each step is propagated **exactly**: the loop is linear with constant coefficients, so
    /// the state `x = (dv, dth)` evolves as `x(t + h) = Phi x(t) + A^-1 (Phi - I) u` with
    /// `A = [[0, -g], [1/R, 0]]`, constant forcing `u = (b_a, eps)` and `Phi = exp(A h)`:
    ///
    /// ```text
    /// Phi = [[cos wh,            -(g/w) sin wh],
    ///        [sin wh / (R w),     cos wh      ]],   A^-1 = [[0, R], [-1/g, 0]]
    /// ```
    ///
    /// The white noise is added once per step with variance `PSD * h`, which is accurate
    /// to first order in `w h` (about 1e-3 for 1 s steps).
    pub fn velocity_errors<R: Rng + ?Sized>(
        &self,
        biases: InsBiases,
        dt: f64,
        n_steps: usize,
        with_noise: bool,
        rng: &mut R,
    ) -> Vec<Vec2> {
        let (g, r, w) = (GRAVITY_M_S2, EARTH_RADIUS_M, schuler_omega());
        let (c, s) = ((w * dt).cos(), (w * dt).sin());
        let phi = [[c, -(g / w) * s], [s / (r * w), c]];
        let vrw = self.spec.vrw_m_per_s_per_sqrt_s() * dt.sqrt();
        let arw = self.spec.arw_rad_per_sqrt_s() * dt.sqrt();

        // One axis at a time; the closure captures the shared matrices by reference.
        let axis = |b_a: f64, eps: f64, rng: &mut R| -> Vec<f64> {
            // Particular part A^-1 (Phi - I) u, the same every step.
            let du = [
                (phi[0][0] - 1.0) * b_a + phi[0][1] * eps,
                phi[1][0] * b_a + (phi[1][1] - 1.0) * eps,
            ];
            let forced = [r * du[1], -du[0] / g];
            let (mut dv, mut dth) = (0.0, 0.0);
            (0..n_steps)
                .map(|_| {
                    let next_dv = phi[0][0] * dv + phi[0][1] * dth + forced[0];
                    let next_dth = phi[1][0] * dv + phi[1][1] * dth + forced[1];
                    dv = next_dv;
                    dth = next_dth;
                    if with_noise {
                        dv += vrw * rng.sample::<f64, _>(StandardNormal);
                        dth += arw * rng.sample::<f64, _>(StandardNormal);
                    }
                    dv
                })
                .collect()
        };
        let ex = axis(biases.accel_m_s2.x, biases.gyro_rad_s.x, rng);
        let ey = axis(biases.accel_m_s2.y, biases.gyro_rad_s.y, rng);
        ex.into_iter()
            .zip(ey)
            .map(|(x, y)| Vec2::new(x, y))
            .collect()
    }

    /// Approximate 1-sigma velocity error per axis at time `t` (m/s), for sizing a filter.
    ///
    /// Bias terms are exact: `sigma_a sin(wt) / w` and `R sigma_eps (1 - cos wt)`. Noise terms
    /// use the short-time forms `VRW sqrt(t)` and `g ARW sqrt(t^3 / 3)` (tilt random walk,
    /// integrated once through gravity). Terms are independent, so they add in quadrature.
    pub fn velocity_error_sigma(&self, t: f64) -> f64 {
        let w = schuler_omega();
        let accel = self.spec.accel_bias_m_s2() * (w * t).sin() / w;
        let gyro = EARTH_RADIUS_M * self.spec.gyro_bias_rad_s() * (1.0 - (w * t).cos());
        let vrw = self.spec.vrw_m_per_s_per_sqrt_s() * t.sqrt();
        let arw = GRAVITY_M_S2 * self.spec.arw_rad_per_sqrt_s() * (t.powi(3) / 3.0).sqrt();
        (accel * accel + gyro * gyro + vrw * vrw + arw * arw).sqrt()
    }
}

/// The INS error model used by a simulated flight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InsModel {
    /// Constant velocity bias (synthetic test model).
    ConstantBias(ConstantBiasIns),
    /// Unaided INS on a real IMU, with Schuler-loop error growth.
    Schuler(SchulerIns),
}

/// `From` lets callers write `ConstantBiasIns { .. }.into()` where an `InsModel` is expected.
impl From<ConstantBiasIns> for InsModel {
    fn from(ins: ConstantBiasIns) -> Self {
        InsModel::ConstantBias(ins)
    }
}

impl InsModel {
    /// INS velocity errors at `t = dt, ..., n_steps dt` for one flight. The constant-bias
    /// model uses no random numbers, so it leaves `rng` untouched.
    pub fn velocity_errors<R: Rng + ?Sized>(
        &self,
        dt: f64,
        n_steps: usize,
        rng: &mut R,
    ) -> Vec<Vec2> {
        match self {
            InsModel::ConstantBias(ins) => vec![ins.velocity_bias; n_steps],
            InsModel::Schuler(ins) => {
                let biases = ins.draw_biases(rng);
                ins.velocity_errors(biases, dt, n_steps, true, rng)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    const DT: f64 = 1.0;
    const N: usize = 600;

    fn rng() -> Xoshiro256PlusPlus {
        Xoshiro256PlusPlus::seed_from_u64(3)
    }

    fn ins() -> SchulerIns {
        SchulerIns {
            spec: ImuSpec::HG4930_CA51,
        }
    }

    /// Pass: |simulated - closed form| <= 1e-3 * max |closed form| at every step.
    fn assert_matches(sim: &[Vec2], exact: impl Fn(f64) -> Vec2) {
        let scale = (1..=N)
            .map(|k| exact(k as f64 * DT).norm())
            .fold(0.0, f64::max);
        for (k, v) in sim.iter().enumerate() {
            let e = exact((k + 1) as f64 * DT);
            assert!((*v - e).norm() <= 1e-3 * scale, "step {k}: {v:?} vs {e:?}");
        }
    }

    #[test]
    fn accelerometer_bias_alone_matches_closed_form() {
        // Closed form (check by differentiating): dv = (b/w) sin wt, dth = (b/g)(1 - cos wt).
        let b = Vec2::new(3e-3, -1e-3);
        let biases = InsBiases {
            accel_m_s2: b,
            gyro_rad_s: Vec2::default(),
        };
        let sim = ins().velocity_errors(biases, DT, N, false, &mut rng());
        let w = schuler_omega();
        assert_matches(&sim, |t| b * ((w * t).sin() / w));
    }

    #[test]
    fn gyro_bias_alone_matches_closed_form() {
        // Closed form: dth = (eps/w) sin wt, dv = -R eps (1 - cos wt).
        let eps = Vec2::new(5e-6, 2e-6);
        let biases = InsBiases {
            accel_m_s2: Vec2::default(),
            gyro_rad_s: eps,
        };
        let sim = ins().velocity_errors(biases, DT, N, false, &mut rng());
        let w = schuler_omega();
        assert_matches(&sim, |t| eps * (-EARTH_RADIUS_M * (1.0 - (w * t).cos())));
    }

    #[test]
    fn monte_carlo_spread_matches_sizing_formula() {
        // Pass: std of the simulated velocity error over 2000 flights at t = 600 s, per axis,
        // is within 5% of `velocity_error_sigma(600)` (standard error of a std estimate with
        // 2000 samples is about 1.6%).
        let ins = ins();
        let mut rng = rng();
        let runs = 2000;
        let (mut sx, mut sy) = (0.0, 0.0);
        for _ in 0..runs {
            let biases = ins.draw_biases(&mut rng);
            let last = *ins
                .velocity_errors(biases, DT, N, true, &mut rng)
                .last()
                .unwrap();
            sx += last.x * last.x;
            sy += last.y * last.y;
        }
        let expected = ins.velocity_error_sigma(N as f64 * DT);
        for s in [sx, sy] {
            let std = (s / runs as f64).sqrt();
            assert!(
                (std / expected - 1.0).abs() < 0.05,
                "std {std:.3} vs {expected:.3} m/s"
            );
        }
    }

    #[test]
    fn constant_bias_model_uses_no_randomness() {
        let model: InsModel = ConstantBiasIns {
            velocity_bias: Vec2::new(0.0, 2.0),
        }
        .into();
        let mut a = rng();
        let errors = model.velocity_errors(DT, 3, &mut a);
        assert_eq!(errors, vec![Vec2::new(0.0, 2.0); 3]);
        assert_eq!(a.random::<u64>(), rng().random::<u64>());
    }
}
