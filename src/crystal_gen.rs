use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use rayon::prelude::*;
use std::collections::BTreeMap;

use libcrystal::{Crystal as crystal_, CrystalGenerator as crystal_gen, Float};

#[pyclass(module = "crystallus")]
#[text_signature = "(spacegroup_num, estimated_volume, estimated_variance, *, min_distance_tolerance, angle_range, angle_tolerance, max_recurrent, n_jobs)"]
pub struct CrystalGenerator {
    _crystal_gen: crystal_gen<'static>,
    _n_jobs: i16,
}

#[pymethods]
impl CrystalGenerator {
    #[new]
    #[args(
        "*",
        min_distance_tolerance = "0.15",
        angle_range = "(30., 150.)",
        angle_tolerance = "20.",
        max_recurrent = "5_000",
        n_jobs = "-1",
        verbose = true
    )]
    fn new(
        spacegroup_num: usize,
        estimated_volume: Float,
        estimated_variance: Float,
        min_distance_tolerance: Float,
        angle_range: (Float, Float),
        angle_tolerance: Float,
        max_recurrent: u16,
        n_jobs: i16,
        verbose: bool,
    ) -> PyResult<Self> {
        let _crystal_gen = crystal_gen::from_spacegroup_num(
            spacegroup_num,
            estimated_volume,
            estimated_variance,
            Some(min_distance_tolerance),
            Some(angle_range),
            Some(angle_tolerance),
            Some(max_recurrent),
            Some(verbose),
        );
        match _crystal_gen {
            Err(e) => Err(exceptions::ValueError::py_err(e.to_string())),
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

    #[getter(min_distance_tolerance)]
    fn min_distance_tolerance(&self) -> PyResult<Float> {
        Ok(self._crystal_gen.min_distance_tolerance)
    }

    #[getter(max_recurrent)]
    fn max_recurrent(&self) -> PyResult<u16> {
        Ok(self._crystal_gen.max_recurrent)
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
    #[text_signature = "($self, check_distance, **cfg)"]
    #[args(cfg = "**")]
    fn gen_one(
        &self,
        py: Python<'_>,
        check_distance: bool,
        cfg: Option<&PyDict>,
    ) -> PyResult<PyObject> {
        match cfg {
            Some(cfg) => {
                let mut cfg: BTreeMap<&str, Vec<&str>> = cfg.extract()?;
                let mut elements: Vec<&str> = Vec::new();
                let mut wyckoff_letters: Vec<&str> = Vec::new();
                for (elem, letter) in cfg.iter_mut() {
                    elements.append(&mut vec![elem; letter.len()]);
                    wyckoff_letters.append(letter);
                }
                let cry = self
                    ._crystal_gen
                    .gen(&elements, &wyckoff_letters, Some(check_distance));

                match cry {
                    Err(e) => Err(exceptions::ValueError::py_err(e.to_string())),
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
                return Err(exceptions::ValueError::py_err(
                    "no configurations for generation",
                ));
            }
        }
    }

    #[args(cfgs = "*")]
    #[text_signature = "($self, size, check_distance, /, *cfgs)"]
    fn gen_many(
        &self,
        py: Python<'_>,
        size: i32,
        check_distance: bool,
        cfgs: &PyTuple,
    ) -> PyResult<PyObject> {
        let mut cfgs: Vec<BTreeMap<&str, Vec<&str>>> =
            match cfgs.extract() {
                Ok(m) => m,
                Err(_) => return Err(exceptions::ValueError::py_err(
                    "can not converting `cfg` into dict, make sure the `cfgs` are tuple of dicts",
                )),
            };
        // parallel using rayon
        if self._n_jobs > 0 {
            std::env::set_var("RAYON_NUM_THREADS", self._n_jobs.to_string());
        }

        let mut ret: Vec<crystal_> = Vec::new();
        match cfgs.len() {
            0 => {
                return Err(exceptions::ValueError::py_err(
                    "no configurations for generation",
                ));
            }
            1 => {
                let mut elements: Vec<&str> = Vec::new();
                let mut wyckoff_letters: Vec<&str> = Vec::new();
                for (elem, letter) in cfgs[0].iter_mut() {
                    elements.append(&mut vec![elem; letter.len()]);
                    wyckoff_letters.append(letter);
                }

                //Do works
                ret.append(&mut py.allow_threads(move || {
                    (0..size)
                        .into_par_iter()
                        .map(|_| {
                            self._crystal_gen
                                .gen(&elements, &wyckoff_letters, Some(check_distance))
                        })
                        .filter_map(Result::ok)
                        .collect::<Vec<crystal_>>()
                }));
            }
            _ => {
                for cfg in cfgs.iter_mut() {
                    let mut elements: Vec<&str> = Vec::new();
                    let mut wyckoff_letters: Vec<&str> = Vec::new();
                    for (elem, letter) in cfg.iter_mut() {
                        elements.append(&mut vec![elem; letter.len()]);
                        wyckoff_letters.append(letter);
                    }

                    //Do works
                    ret.append(&mut py.allow_threads(move || {
                        (0..size)
                            .into_par_iter()
                            .map(|_| {
                                self._crystal_gen.gen(
                                    &elements,
                                    &wyckoff_letters,
                                    Some(check_distance),
                                )
                            })
                            .filter_map(Result::ok)
                            .collect::<Vec<crystal_>>()
                    }));
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
