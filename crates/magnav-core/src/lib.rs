//! # magnav-core
//!
//! Core library for simulating magnetic anomaly navigation (MagNav) over the UK.
//!
//! ## Conventions (used everywhere in this crate)
//! - Positions are planar coordinates in **metres**: `x` = east, `y` = north.
//!   Real maps will use the British National Grid (EPSG:27700).
//! - Time is in **seconds**, velocities in **m/s**.
//! - Magnetic field values are total-field anomalies in **nanotesla (nT)**.

pub mod map;

pub use map::{FieldMap, MapError, MapGrid, MapMetadata};
