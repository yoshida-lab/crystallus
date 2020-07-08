// Copyright 2020 TsumiNa. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#[macro_use]
extern crate lazy_static;

mod crystal_gen;
mod wyckoff_cfg;
mod wyckoff_pos;

pub use crystal_gen::*;
pub use wyckoff_cfg::*;
pub use wyckoff_pos::*;

use std::collections::HashMap;

#[cfg(feature = "f32")]
pub type Float = f32;
#[cfg(not(feature = "f32"))]
pub type Float = f64;

/// Space group type for all 230 space groups.
pub const SPG_TYPES: [char; 230] = [
    'P', 'P', 'P', 'P', 'C', 'P', 'P', 'C', 'C', 'P', 'P', 'C', 'P', 'P', 'C', 'P', 'P', 'P', 'P',
    'C', 'C', 'F', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'C', 'C', 'C', 'A',
    'A', 'A', 'A', 'F', 'F', 'I', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P',
    'P', 'P', 'P', 'P', 'P', 'C', 'C', 'C', 'C', 'C', 'C', 'F', 'F', 'I', 'I', 'I', 'I', 'P', 'P',
    'P', 'P', 'I', 'I', 'P', 'I', 'P', 'P', 'P', 'P', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P',
    'P', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'I', 'I', 'I', 'I', 'P', 'P', 'P', 'P',
    'P', 'P', 'P', 'P', 'I', 'I', 'I', 'I', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P',
    'P', 'P', 'P', 'P', 'P', 'I', 'I', 'I', 'I', 'P', 'P', 'P', 'R', 'P', 'R', 'P', 'P', 'P', 'P',
    'P', 'P', 'R', 'P', 'P', 'P', 'P', 'R', 'R', 'P', 'P', 'P', 'P', 'R', 'R', 'P', 'P', 'P', 'P',
    'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P',
    'P', 'P', 'P', 'P', 'P', 'F', 'I', 'P', 'I', 'P', 'P', 'F', 'F', 'I', 'P', 'I', 'P', 'P', 'F',
    'F', 'I', 'P', 'P', 'I', 'P', 'F', 'I', 'P', 'F', 'I', 'P', 'P', 'P', 'P', 'F', 'F', 'F', 'F',
    'I', 'I',
];

/// Wyckoff table for all 230 space groups.
const WYCKOFFS: &'static str = std::include_str!("wyckoffs.json");
lazy_static! {
    static ref WY: Vec<HashMap<String, (usize, bool, String)>> =
        serde_json::from_str(WYCKOFFS).unwrap();
}
