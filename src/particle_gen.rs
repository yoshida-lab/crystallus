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

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use pyo3::types::PyTuple;
use rayon::prelude::*;

use libcrystal::{Coord, Float, WyckoffPos as wy_pos};

#[pyclass(module = "crystallus")]
pub struct ParticleGenerator {
    _wy_pos: wy_pos,
}

#[pymethods]
impl ParticleGenerator {
    #[new]
    #[args(shifts = "*")]
    fn new(s: &str, shifts: &PyTuple) -> PyResult<Self> {
        match shifts.len() {
            l if l > 0 => {
                let mut ret: Vec<Coord> = Vec::new();
                for shift in shifts.iter() {
                    let shift: Vec<Float> = shift.extract()?;
                    if shift.len() != 3 {
                        return Err(PyValueError::new_err(format!(
                            "len of {:?} must be 3",
                            shift
                        )));
                    }
                    ret.push([shift[0], shift[1], shift[2]]);
                }
                let _wy_pos = wy_pos::from_str_and_shifts(s, ret);
                match _wy_pos {
                    Err(e) => Err(PyValueError::new_err(e.to_string())),
                    Ok(w) => return Ok(ParticleGenerator { _wy_pos: w }),
                }
            }
            _ => {
                let _wy_pos = s.parse::<wy_pos>();
                match _wy_pos {
                    Err(e) => Err(PyValueError::new_err(e.to_string())),
                    Ok(w) => return Ok(ParticleGenerator { _wy_pos: w }),
                }
            }
        }
    }

    fn random_gen(&self, py: Python<'_>) -> PyResult<PyObject> {
        return Ok(self._wy_pos.random_gen().into_raw_vec().into_py(py));
    }

    #[args(only_return_len = "false")]
    fn gen_many(&self, py: Python<'_>, size: i32, only_return_len: bool) -> PyResult<PyObject> {
        let ret = (0..size)
            .into_iter()
            .map(|_| self._wy_pos.random_gen().into_raw_vec())
            .collect::<Vec<Vec<Float>>>();

        if only_return_len {
            return Ok(ret.len().into_py(py));
        } else {
            return Ok(ret.into_py(py));
        }
    }

    #[args(only_return_len = "false")]
    fn par_gen_many(&self, py: Python<'_>, size: i32, only_return_len: bool) -> PyResult<PyObject> {
        let ret = py.allow_threads(move || {
            (0..size)
                .into_par_iter()
                .map(|_| self._wy_pos.random_gen().into_raw_vec())
                .collect::<Vec<Vec<Float>>>()
        });

        if only_return_len {
            return Ok(ret.len().into_py(py));
        } else {
            return Ok(ret.into_py(py));
        }
    }
}
