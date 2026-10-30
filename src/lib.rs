//! Rust numerical implementations for QRK-Sylv and its four named baselines.
//!
//! Published quantum baselines are represented by their mathematical solution-map
//! emulators, NOT by claimed gate-level implementations of the original circuits.
pub mod complex;
pub mod matrix;
pub mod rng;
pub mod sampling;
pub mod linalg;
pub mod solvers;
pub mod fixtures;
pub mod clock;
