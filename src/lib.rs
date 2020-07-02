use pyo3::exceptions;
use pyo3::prelude::*;
use pyo3::types::PyList;

mod crystal_gen;
mod particle_gen;
mod wyckoff_cfg_gen;

use crate::crystal_gen::CrystalGenerator;
use crate::particle_gen::ParticleGenerator;
use crate::wyckoff_cfg_gen::WyckoffCfgGenerator;
use libcrystal::{lll_reduce as _lll, pbc_all_distances as _pbc, Float};

#[pymodule(crystallus)]
fn crystallus(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<ParticleGenerator>()?;
    m.add_class::<CrystalGenerator>()?;
    m.add_class::<WyckoffCfgGenerator>()?;

    #[pyfn(m, "lll_reduce", "*", delta = "0.75")]
    #[text_signature = "(basis, delta)"]
    fn lll_reduce(_py: Python, basis: &PyList, delta: Float) -> PyResult<(Vec<Float>, Vec<Float>)> {
        let basis: Vec<Vec<Float>> = basis.extract()?;
        let basis: Vec<[Float; 3]> = basis.iter().map(|x| [x[0], x[1], x[2]]).collect();
        let (basis, mapping) = _lll(&basis, delta);
        Ok((basis, mapping))
    }

    #[pyfn(m, "pbc_all_distances")]
    #[text_signature = "(lattice, frac_coords)"]
    fn pbc_all_distances(
        _py: Python,
        lattice: &PyList,
        frac_coords: &PyList,
    ) -> PyResult<Vec<Vec<Float>>> {
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
            Err(e) => Err(exceptions::ValueError::py_err(format!("{}", e))),
        }
    }

    Ok(())
}
