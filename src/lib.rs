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

use pyo3::{prelude::*, py_run};

mod crystal_gen;
mod particle_gen;
mod utils;
mod wyckoff_cfg_gen;

use crate::crystal_gen::CrystalGenerator;
use crate::particle_gen::ParticleGenerator;
use crate::wyckoff_cfg_gen::WyckoffCfgGenerator;

#[pymodule]
fn core(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();

    // register classes
    m.add_class::<ParticleGenerator>()?;
    m.add_class::<CrystalGenerator>()?;
    m.add_class::<WyckoffCfgGenerator>()?;

    // register functions
    let utils_mod = PyModule::new(py, "utils")?;
    utils::register(py, utils_mod)?;

    // Note that this does not define a package, so this won’t allow Python code to
    // directly import submodules by using from parent_module import child_module.
    // For more information, see https://pyo3.rs/latest/module.html#python-submodules
    py_run!(
        py,
        utils_mod,
        "import sys; sys.modules['crystallus.core.utils'] = utils_mod"
    );

    // add submodules
    m.add_submodule(utils_mod)?;
    Ok(())
}
