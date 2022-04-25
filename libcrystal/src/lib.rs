// Copyright 2021 TsumiNa
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#[macro_use]
extern crate lazy_static;

mod generator;
mod structure;
pub mod utils;
mod wyckoff_cfg;
mod wyckoff_pos;

pub use generator::*;
pub use structure::*;
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
const WYCKOFFS: &'static str = std::include_str!("external/wyckoffs.json");
lazy_static! {
    static ref WY: Vec<HashMap<String, (usize, bool, String)>> =
        serde_json::from_str(WYCKOFFS).unwrap();
}

// Covalent radius for element H (Z=1) to Cm (Z=96)
const COVALENT_RADIUS: &'static str = std::include_str!("external/covalent_radius.json");
lazy_static! {
    static ref RADIUS: HashMap<String, Float> = serde_json::from_str(COVALENT_RADIUS).unwrap();
}
