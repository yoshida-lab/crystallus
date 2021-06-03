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

use crate::Float;
use ndarray::Array2;
use std::{error, fmt};

#[derive(Clone, Debug)]
pub struct RandomGeneratorOption {
    pub angle_range: (Float, Float),
    pub angle_tolerance: Float,
    pub max_attempts_number: u16,
    pub verbose: bool,
}

impl Default for RandomGeneratorOption {
    fn default() -> Self {
        Self {
            angle_range: (30., 150.),
            angle_tolerance: 20.,
            max_attempts_number: 5_000,
            verbose: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TemplateBaseGeneratorOption {
    pub angle_range: (Float, Float),
    pub angle_tolerance: Float,
    pub lattice: Vec<Float>,
    pub empirical_coords: Vec<(String, Vec<Float>)>,
    pub empirical_coords_variance: Float,
    pub empirical_coords_sampling_rate: Float,
    pub empirical_coords_loose_sampling: bool,
    pub max_attempts_number: u16,
    pub verbose: bool,
}

impl Default for TemplateBaseGeneratorOption {
    fn default() -> Self {
        Self {
            angle_range: (30., 150.),
            angle_tolerance: 20.,
            lattice: vec![0.; 9],
            empirical_coords: Vec::new(),
            empirical_coords_variance: 0.01,
            empirical_coords_sampling_rate: 1.,
            empirical_coords_loose_sampling: true,
            max_attempts_number: 5_000,
            verbose: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CrystalGeneratorError(pub String);

impl fmt::Display for CrystalGeneratorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "crystal generator error: `{}`", self.0)
    }
}

impl error::Error for CrystalGeneratorError {}

pub type LatticeFn =
    Box<dyn Fn() -> Result<(Array2<Float>, Float), CrystalGeneratorError> + Send + Sync>;

#[derive(Debug, Clone, PartialEq)]
pub struct Crystal {
    pub spacegroup_num: usize,
    pub volume: Float,
    pub lattice: Array2<Float>,
    pub particles: Array2<Float>,
    pub elements: Vec<String>,
    pub wyckoff_letters: Vec<String>,
}
