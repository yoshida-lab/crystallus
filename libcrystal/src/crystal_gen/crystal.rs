// Copyright 2020 TsumiNa. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::Float;
use ndarray::Array2;

#[derive(Debug, Clone, PartialEq)]
pub struct Crystal {
    pub spacegroup_num: usize,
    pub volume: Float,
    pub lattice: Array2<Float>,
    pub particles: Array2<Float>,
    pub elements: Vec<String>,
    pub wyckoff_letters: Vec<String>,
}
