//! # magnav-core
//!
//! Core library for simulating magnetic anomaly navigation (MagNav) over the UK.
//!
//! ## Conventions (used everywhere in this crate)
//! - Positions are planar coordinates in **metres**: `x` = east, `y` = north.
//!   Real maps will use the British National Grid (EPSG:27700).
//! - Time is in **seconds**, velocities in **m/s**.
//! - Magnetic field values are total-field anomalies in **nanotesla (nT)**.

pub mod filter;
pub mod map;
pub mod sensor;
pub mod sim;
pub mod synthetic;
pub mod trajectory;
pub mod vec2;

pub use filter::{ExecMode, FilterConfig, FilterError, ParticleFilter, StepSummary};
pub use map::{Bounds, FieldMap, MapError, MapGrid, MapMetadata};
pub use sim::{Scenario, SimError, StepRecord, run_scenario};
pub use vec2::Vec2;
