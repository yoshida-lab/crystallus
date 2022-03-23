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

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PySequence;

use libcrystal::{
    wrap::{lll_reduce as _lll, pbc_all_distances as _pbc},
    Float,
};

// register functions
pub(crate) fn register(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(lll_reduce, m)?)?;
    m.add_function(wrap_pyfunction!(pbc_all_distances, m)?)?;
    Ok(())
}

#[pyfunction("*", delta = "0.75")]
#[pyo3(name = "lll_reduce")]
#[pyo3(text_signature = "(basis, delta)")]
fn lll_reduce(basis: &PySequence, delta: Float) -> PyResult<(Vec<Float>, Vec<Float>)> {
    let basis: Vec<Vec<Float>> = basis.extract()?;
    let basis: Vec<[Float; 3]> = basis.iter().map(|x| [x[0], x[1], x[2]]).collect();
    let (basis, mapping) = _lll(&basis, delta);
    Ok((basis, mapping))
}

#[pyfunction("*", delta = "0.75")]
#[pyo3(name = "pbc_all_distances")]
#[pyo3(text_signature = "(lattice, frac_coords)")]
fn pbc_all_distances(lattice: &PySequence, frac_coords: &PySequence) -> PyResult<Vec<Vec<Float>>> {
    let lattice: Vec<Vec<Float>> = lattice.extract()?;
    let frac_coords: Vec<Vec<Float>> = frac_coords.extract()?;
    let lattice: Vec<[Float; 3]> = lattice.iter().map(|x| [x[0], x[1], x[2]]).collect();
    let frac_coords: Vec<[Float; 3]> = frac_coords.iter().map(|x| [x[0], x[1], x[2]]).collect();
    let ret = _pbc(&lattice, &frac_coords);
    match ret {
        Ok(d) => {
            let chunk_size = (d.len() as Float).sqrt() as usize;
            let mut ret_ = Vec::new();
            for chunk in d.chunks(chunk_size) {
                ret_.push(chunk.to_vec());
            }
            Ok(ret_)
        }
        Err(e) => Err(PyValueError::new_err(format!("{}", e))),
    }
}
