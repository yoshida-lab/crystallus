// Copyright 2020 TsumiNa. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Float;
use ndarray::Array2;

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
#[derive(Debug, Clone, PartialEq)]
pub struct Crystal {
    pub spacegroup_num: usize,
    pub volume: Float,
    pub lattice: Array2<Float>,
    pub particles: Array2<Float>,
    pub elements: Vec<String>,
    pub wyckoff_letters: Vec<String>,
}
