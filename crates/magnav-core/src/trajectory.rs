//! Vehicle truth trajectories and the inertial navigation system (INS) error model.
//!
//! Phase 1 uses the simplest useful case: straight, constant-velocity flight, and an INS whose
//! velocity has a constant bias, so its dead-reckoned position error grows linearly in time.
//! Realistic INS error models (bias instability, random walk) come in Phase 2, with
//! datasheet-sourced parameters.

use crate::vec2::Vec2;

/// Straight-line flight at constant velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StraightLine {
    /// Position at t = 0 (m).
    pub start: Vec2,
    /// Constant velocity (m/s).
    pub velocity: Vec2,
}

impl StraightLine {
    /// True position at time `t` (s).
    pub fn position_at(&self, t: f64) -> Vec2 {
        self.start + self.velocity * t
    }
}

/// An INS that reports `true velocity + constant bias`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ConstantBiasIns {
    /// Velocity bias (m/s).
    pub velocity_bias: Vec2,
}

impl ConstantBiasIns {
    /// Velocity the INS reports when the vehicle truly moves at `true_velocity`.
    pub fn reported_velocity(&self, true_velocity: Vec2) -> Vec2 {
        true_velocity + self.velocity_bias
    }

    /// Dead-reckoned position after `t` seconds of straight flight, starting from `ins_start`
    /// (the INS's own belief about the start, which may itself be wrong).
    pub fn dead_reckon(&self, ins_start: Vec2, truth: &StraightLine, t: f64) -> Vec2 {
        ins_start + self.reported_velocity(truth.velocity) * t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn straight_line_moves_at_constant_velocity() {
        let path = StraightLine {
            start: Vec2::new(1000.0, 2000.0),
            velocity: Vec2::new(100.0, -50.0),
        };
        assert_eq!(path.position_at(0.0), path.start);
        assert_eq!(path.position_at(10.0), Vec2::new(2000.0, 1500.0));
    }

    #[test]
    fn ins_error_grows_linearly_with_bias() {
        let truth = StraightLine {
            start: Vec2::new(0.0, 0.0),
            velocity: Vec2::new(100.0, 0.0),
        };
        let ins = ConstantBiasIns {
            velocity_bias: Vec2::new(0.0, 2.0),
        };
        for t in [0.0, 60.0, 600.0] {
            let error = ins.dead_reckon(truth.start, &truth, t) - truth.position_at(t);
            assert_eq!(error, Vec2::new(0.0, 2.0 * t));
        }
    }
}
