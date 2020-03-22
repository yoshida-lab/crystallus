# Crystallus

[![Actions Status](https://github.com/yoshida-lab/crystallus/workflows/tests/badge.svg)](https://github.com/yoshida-lab/crystallus/actions)

A Rust library for crystal structure generation. The python bindings are provided using [pyo3](https://github.com/PyO3/pyo3).

## Installation

To install this package, you have to prepare an environment which has python development toolchain and Rust compiler toolchain.
We support Python 3.6, 3.7, and 3.8 on macOS and Linux (**windows is not supported**).

1. install dev-python

   We recommend you to use conda to build your python environment. We have provided you with a preset conda env file to smooth your installation.
   Suppose we will use Python 3.7 and Rust 1.42.0-nightly 2020-03-06. Also, we suppose you have been installed conda and be familiar with conda commands.
   The following commands will check your conda installation, then create a new environment named `crystallus`, and install packages which are listed in our preset environment fils.

   ```bash
   $> conda -V  # will return your conda version if conda installation is ok
   $> conda create -n crystallus python=3.7  # create a new environment with python3.7 and name it *crystallus*.
   $> conda env update -n crystallus -f env/environment.yml  # install packages which are listed in `environment.yml` file.
   ```

   The environment's name can be anything you liked but without space in the name.

2. install rust toolchain

   As we said, we use [pyo3](https://github.com/PyO3/pyo3) to provide you the python bindings.
   Based on the official document of pyo3, you have to use nightly build Rust. The minimum required Rust version is 1.42.0-nightly 2020-03-06.
   First, you have to install the [rustup](https://www.rust-lang.org/tools/install) tools into your system. This tool enables you to manage your Rust installation locally.
   If you have done, use the following commands to set up your rust compiler toolchain.

   ```bash
    $> rustup toolchain add nightly-2020-03-06
    $> rustup default nightly-2020-03-06
   ```

3. compile and install crystallus

   This is the last step. We will use [maturin](https://github.com/PyO3/maturin) to build our rust codes into the native dynamic library, and wrap it into a python module.

   ```bash
    $> conda activate crystallus  # activate your python environment by name
    $> pip install -U maturin  # install maturin using pip
    $> cd <path/to/crystallus>  # chang to the root of crystallus codes
    $> maturin build -i $(which python) --release --no-sdist --strip  # build package
   ```

   If everything goes well, you will find a `crystallus-*.whl` file under `target/wheels` directory.
   Then you can use pip to install it into your python environment (in this example, it is the newly created environment named `crystallus`).

   ```bash
   $> pip install -U target/wheels/crystallus-*.whl
   ```

## Usage

Please see examples.
