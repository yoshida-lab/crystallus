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

#[derive(Clone, Debug)]
pub struct BaseGeneratorOption {
    pub variance_of_volume: Float,
    pub angle_range: (Float, Float),
    pub angle_tolerance: Float,
    pub max_attempts_number: u16,
    pub verbose: bool,
}

impl Default for BaseGeneratorOption {
    fn default() -> Self {
        Self {
            variance_of_volume: 10.,
            angle_range: (30., 150.),
            angle_tolerance: 20.,
            max_attempts_number: 5_000,
            verbose: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct EmpiricalGeneratorOption {
    pub base_option: BaseGeneratorOption,
    pub lattice: Vec<Float>,
    pub empirical_coords: Vec<(String, Vec<Float>)>,
    pub empirical_coords_variance: Float,
    pub empirical_coords_sampling_rate: Float,
    pub empirical_coords_loose_sampling: bool,
}

impl Default for EmpiricalGeneratorOption {
    fn default() -> Self {
        Self {
            base_option: BaseGeneratorOption::default(),
            lattice: vec![0.; 9],
            empirical_coords: Vec::new(),
            empirical_coords_variance: 0.01,
            empirical_coords_sampling_rate: 1.,
            empirical_coords_loose_sampling: true,
        }
    }
}
