# Crystallus

[![Actions Status](https://github.com/yoshida-lab/crystallus/workflows/tests/badge.svg)](https://github.com/yoshida-lab/crystallus/actions)

A Rust library for crystal structure generation.
This package is heavily developed on rust to guarantee the performance and provide a python interface to easy the use.

## Installation

To install this package, you have to prepare an environment which has python development toolchain and Rust compiler toolchain.
We support Python 3.6, 3.7, and 3.8 on macOS and Linux (**windows is not supported**).
We highly recommend you to use [miniconda](https://docs.conda.io/en/latest/miniconda.html) and our preset conda env file (`env/environment.yml`) to build your python environment.

1. install dev-python

    Suppose we will use Python 3.7 and Rust 1.46.0-nightly 2020-06-26.
    The following commands will check your conda installation. If everything is ok, then create a new environment named `crystallus` and install packages using `env/environment.yml` environment file.

    ```bash
    $> conda -V  # will return your conda version if conda installation is ok
    $> conda create -n crystallus python=3.7  # create a new environment with python3.7 and name it *crystallus*.
    $> conda env update -n crystallus -f env/environment.yml  # install packages which are listed in `environment.yml` file.
    ```

    The environment's name can be anything you liked without space in the name string.

2. install rust toolchain

    We use [pyo3](https://github.com/PyO3/pyo3) to provide you the python bindings.
    Based on the official document of pyo3, you have to use nightly build Rust. The tested Rust version is 1.44.0-nightly 2020-06-26.
    We will install the [rustup](https://www.rust-lang.org/tools/install) tools into our system and set the default rust compiler using the following commands.

    ```bash
     $> rustup toolchain add nightly-2020-06-26
     $> rustup default nightly-2020-06-26
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
