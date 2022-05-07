// Copyright 2022 TsumiNa
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
mod base;
mod empirical;
mod options;

pub use self::base::*;
pub use self::empirical::*;
pub use self::options::*;

use crate::Float;
use crate::{Crystal, CrystalGeneratorError};

pub trait Generator<T> {
    /// Create a [`Generator`] instance from `spacegroup_num` and `options`.
    fn from_spacegroup_num(
        spacegroup_num: usize,
        volume_of_cell: Float,
        options: Option<T>,
    ) -> Result<Self, CrystalGeneratorError>
    where
        Self: Sized;

    /// Return crystal structure.
    fn gen(
        &self,
        elements: &Vec<String>,
        wyckoff_letters: &Vec<String>,
        check_distance: Option<bool>,
        distance_scale_factor: Option<Float>,
    ) -> Result<Crystal, CrystalGeneratorError>;
}
