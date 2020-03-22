use pyo3::prelude::*;

mod crystal_gen;
mod particle_gen;
mod wyckoff_cfg_gen;

use crate::crystal_gen::CrystalGenerator;
use crate::particle_gen::ParticleGenerator;
use crate::wyckoff_cfg_gen::WyckoffCfgGenerator;

#[pymodule(crystallus)]
fn crystallus(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<ParticleGenerator>()?;
    m.add_class::<CrystalGenerator>()?;
    m.add_class::<WyckoffCfgGenerator>()?;

    Ok(())
}
