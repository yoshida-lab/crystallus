# Crystallus

A Rust library for crystal generation. The python binding is provided using pyo3.

## Installation

To install this package, you need to prepare the dev-python and Rust compilation toolchain.

1. install dev-python
   We recommend you using conda to build your python environment. The simplest way is using our preset conda env file.
   Use the following commands to check your conda installation is ok, and install our preset environment.

   ```bash
   $> conda -V  # will return your conda version if conda installation is ok
   $> conda create -n crystallus python=3.7  # create a new environment with python3.7 and name it *crystallus*.
   $> conda env update -n crystallus -f env/environment.yml  # install packages which are listed in `environment.yml` file.
   ```

   The environment name can be anything you liked but without space in the name. Python version 3.6, 3.7, and 3.8 are supported.

2. install rust toolchain
   We use [pyo3](https://github.com/PyO3/pyo3) to provide the python bindings.
   Based on the official document of pyo3, you have to use nightly build Rust.
   First, you have to install the [rustup](https://www.rust-lang.org/tools/install) tools into your system.
   If you have done, the `rustup` command is ready for use. Use the following commands to set up your rust toolchain.

   ```bash
    $> rustup toolchain add nightly-2020-03-06
    $> rustup default nightly-2020-03-06
   ```

3. compile and install crystallus
   This is the last step. We will use [maturin](https://github.com/PyO3/maturin) to build our rust codes into a loadable python module.

   ```bash
    $> conda activate crystallus  # activate your python environment by name
    $> pip install -U maturin  # install maturin using pip
    $> cd <path/to/crystallus>  # chang to the root of crystallus codes
    $> maturin build -i $(which python) --release --no-sdist --strip . # build package
   ```

   If everything goes well, you will find a `crystallus-*.whl` file under `target/wheels` directory. Then using pip to install it into your python environment.

   ```bash
   $> pip install target/wheels/crystallus-*.whl
   ```

## Usage

Please see examples.
