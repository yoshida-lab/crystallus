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
use pyo3::types::{PyDict, PyTuple};
use rayon::prelude::*;
use std::collections::BTreeMap;

use libcrystal::{Crystal as crystal_, CrystalGenerator as crystal_gen, Float};

#[pyclass(module = "crystallus")]
#[text_signature = "(spacegroup_num, estimated_volume, estimated_variance, *, min_distance_tolerance, angle_range, angle_tolerance, max_recurrent, n_jobs)"]
pub struct CrystalGenerator {
    _crystal_gen: crystal_gen,
    _n_jobs: i16,
}

#[pymethods]
impl CrystalGenerator {
    #[new]
    #[args(
        "*",
        angle_range = "(30., 150.)",
        angle_tolerance = "20.",
        max_attempts_number = "5_000",
        empirical_coords = "None",
        empirical_coords_variance = "0.01",
        n_jobs = "-1",
        verbose = true
    )]
    fn new(
        spacegroup_num: usize,
        estimated_volume: Float,
        estimated_variance: Float,
        angle_range: (Float, Float),
        angle_tolerance: Float,
        empirical_coords: Option<&PyTuple>,
        empirical_coords_variance: Float,
        max_attempts_number: u16,
        n_jobs: i16,
        verbose: bool,
    ) -> PyResult<Self> {
        // convert Option<T: FromPyObject> -> Option<D>
        // if T.extract() return Err(e), pass this panic to python side
        let empirical_coords = if let Some(t) = empirical_coords {
            let t: Vec<(String, Vec<Float>)> = t.extract()?;
            Some(t)
        } else {
            None
        };
        let _crystal_gen = crystal_gen::from_spacegroup_num(
            spacegroup_num,
            estimated_volume,
            estimated_variance,
            Some(angle_range),
            Some(angle_tolerance),
            empirical_coords,
            Some(empirical_coords_variance),
            Some(max_attempts_number),
            Some(verbose),
        );
        match _crystal_gen {
            Err(e) => Err(PyValueError::new_err(e.to_string())),
            Ok(w) => {
                return Ok(CrystalGenerator {
                    _crystal_gen: w,
                    _n_jobs: n_jobs,
                })
            }
        }
    }

    #[getter(spacegroup_num)]
    fn spacegroup_num(&self) -> PyResult<usize> {
        Ok(self._crystal_gen.spacegroup_num)
    }

    #[getter(max_attempts_number)]
    fn max_attempts_number(&self) -> PyResult<u16> {
        Ok(self._crystal_gen.max_attempts_number)
    }

    #[getter(n_jobs)]
    fn n_jobs(&self) -> PyResult<i16> {
        Ok(self._n_jobs)
    }

    #[setter(n_jobs)]
    fn set_n_jobs(&mut self, n: i16) -> PyResult<()> {
        self._n_jobs = n;
        Ok(())
    }

    #[getter(verbose)]
    fn verbose(&self) -> PyResult<bool> {
        Ok(self._crystal_gen.verbose)
    }

    #[setter(verbose)]
    fn set_verbose(&mut self, verbose: bool) -> PyResult<()> {
        self._crystal_gen.verbose = verbose;
        Ok(())
    }
    #[text_signature = "($self, check_distance, atomic_distance_tolerance, /, **cfg)"]
    #[args(check_distance = true, atomic_distance_tolerance = "0.1", cfg = "**")]
    fn gen_one(
        &self,
        py: Python<'_>,
        check_distance: bool,
        distance_scale_factor: Float,
        cfg: Option<&PyDict>,
    ) -> PyResult<PyObject> {
        match cfg {
            Some(cfg) => {
                let mut cfg: BTreeMap<String, Vec<String>> = cfg.extract()?;
                let mut elements: Vec<String> = Vec::new();
                let mut wyckoff_letters: Vec<String> = Vec::new();
                for (elem, letter) in cfg.iter_mut() {
                    elements.append(&mut vec![(*elem).clone(); letter.len()]);
                    wyckoff_letters.append(letter);
                }
                let cry = self._crystal_gen.gen(
                    &elements,
                    &wyckoff_letters,
                    Some(check_distance),
                    Some(distance_scale_factor),
                );

                match cry {
                    Err(e) => Err(PyValueError::new_err(e.to_string())),
                    Ok(w) => {
                        let dict = PyDict::new(py);
                        dict.set_item("spacegroup_num", w.spacegroup_num)?;
                        dict.set_item("volume", w.volume)?;
                        dict.set_item(
                            "lattice",
                            w.lattice
                                .into_raw_vec()
                                .chunks(3)
                                .collect::<Vec<&[Float]>>(),
                        )?;
                        dict.set_item("species", w.elements)?;
                        dict.set_item("wyckoff_letters", w.wyckoff_letters)?;
                        dict.set_item(
                            "coords",
                            w.particles
                                .into_raw_vec()
                                .chunks(3)
                                .collect::<Vec<&[Float]>>(),
                        )?;

                        Ok(dict.into_py(py))
                    }
                }
            }
            None => {
                return Err(PyValueError::new_err("no configurations for generation"));
            }
        }
    }

    #[text_signature = "($self, expect_size, max_attempts, check_distance, distance_scale_factor, /, *cfgs)"]
    #[args(
        max_attempts = "None",
        check_distance = true,
        distance_scale_factor = "0.1",
        cfgs = "*"
    )]
    fn gen_many(
        &self,
        py: Python<'_>,
        expect_size: usize,
        max_attempts: Option<usize>,
        check_distance: bool,
        distance_scale_factor: Float,
        cfgs: &PyTuple,
    ) -> PyResult<PyObject> {
        let mut cfgs: Vec<BTreeMap<String, Vec<String>>> =
            match cfgs.extract() {
                Ok(m) => m,
                Err(_) => return Err(PyValueError::new_err(
                    "can not converting `cfg` into dict, make sure the `cfgs` are tuple of dicts",
                )),
            };
        // parallel using rayon
        if self._n_jobs > 0 {
            std::env::set_var("RAYON_NUM_THREADS", self._n_jobs.to_string());
        }

        let max_attempts = max_attempts.unwrap_or(expect_size);
        if max_attempts < expect_size {
            return Err(PyValueError::new_err(
                "`max_attempts` can not be smaller than `expect_size`",
            ));
        }
        let mut ret: Vec<crystal_> = Vec::new();
        match cfgs.len() {
            0 => {
                return Ok(PyTuple::new(py, Vec::<PyDict>::new()).into_py(py));
            }
            1 => {
                let mut counter = expect_size;
                let mut elements: Vec<String> = Vec::new();
                let mut wyckoff_letters: Vec<String> = Vec::new();
                for (elem, letter) in cfgs[0].iter_mut() {
                    elements.append(&mut vec![(*elem).clone(); letter.len()]);
                    wyckoff_letters.append(letter);
                }

                while (ret.len() < expect_size) && (counter <= max_attempts) {
                    //Do works
                    ret.append(&mut py.allow_threads(|| {
                        (0..expect_size)
                            .into_par_iter()
                            .map(|_| {
                                self._crystal_gen.gen(
                                    &elements,
                                    &wyckoff_letters,
                                    Some(check_distance),
                                    Some(distance_scale_factor),
                                )
                            })
                            .filter_map(Result::ok)
                            .collect::<Vec<crystal_>>()
                    }));
                    counter += expect_size;
                }
            }
            _ => {
                let mut ref_ = 0;
                for cfg in cfgs.iter_mut() {
                    let mut counter = expect_size;
                    let mut elements: Vec<String> = Vec::new();
                    let mut wyckoff_letters: Vec<String> = Vec::new();
                    for (elem, letter) in cfg.iter_mut() {
                        elements.append(&mut vec![(*elem).clone(); letter.len()]);
                        wyckoff_letters.append(letter);
                    }

                    while (ret.len() - ref_ < expect_size) && (counter <= max_attempts) {
                        //Do works
                        ret.append(&mut py.allow_threads(|| {
                            (0..expect_size)
                                .into_par_iter()
                                .map(|_| {
                                    self._crystal_gen.gen(
                                        &elements,
                                        &wyckoff_letters,
                                        Some(check_distance),
                                        Some(distance_scale_factor),
                                    )
                                })
                                .filter_map(Result::ok)
                                .collect::<Vec<crystal_>>()
                        }));
                        counter += expect_size;
                    }
                    ref_ = ret.len();
                }
            }
        }

        std::env::set_var("RAYON_NUM_THREADS", "");

        let mut ret_: Vec<PyObject> = Vec::new();
        for crystal in ret {
            let dict = PyDict::new(py);
            dict.set_item("spacegroup_num", crystal.spacegroup_num)?;
            dict.set_item("volume", crystal.volume)?;
            dict.set_item(
                "lattice",
                crystal
                    .lattice
                    .into_raw_vec()
                    .chunks(3)
                    .collect::<Vec<&[Float]>>(),
            )?;
            dict.set_item("species", crystal.elements)?;
            dict.set_item("wyckoff_letters", crystal.wyckoff_letters)?;
            dict.set_item(
                "coords",
                crystal
                    .particles
                    .into_raw_vec()
                    .chunks(3)
                    .collect::<Vec<&[Float]>>(),
            )?;

            ret_.push(dict.into_py(py));
        }
        Ok(PyTuple::new(py, ret_).into_py(py))
    }
}
