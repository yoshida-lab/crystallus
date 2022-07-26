# Crystallus

[![Actions Status](https://github.com/yoshida-lab/crystallus/workflows/tests/badge.svg)](https://github.com/yoshida-lab/crystallus/actions)

A Rust library for crystal structure generation.
This package is heavily developed on rust to guarantee the performance and provide a python interface to easy the use.

## Installation

To install this package, you have to prepare an environment which has python development toolchain and Rust compiler toolchain.
We support Python 3.7, 3.8, and 3.9 on macOS and Linux (**windows is not tested, you should do it on your risks**).
We highly recommend you to use [miniconda](https://docs.conda.io/en/latest/miniconda.html) and our preset conda env file (`env/environment.yml`) to build your python environment.

Suppose we will use the Python 3.8 and the latest stable Rust on Linux.

1. install libopenblas-dev and gfortran onto your Linux system

    We have to compile this package with locally installed lapack packages until now. The final goal is to compile this package in conda.

    ```bash
    $> sudo apt update
    $> sudo apt install -y libopenblas-dev gfortran
    ```

2. install dev-python

    The following commands will check your conda installation. If everything is ok, then create a new environment named `crystallus` and install packages using `env/environment.yml` environment file.

    ```bash
    $> conda -V  # will return your conda version if conda installation is ok
    $> conda create -n crystallus python=3.8  # create a new environment with python3.7 and name it *crystallus*.
    $> conda env update -n crystallus -f env/environment.yml  # install packages which are listed in `environment.yml` file.
    ```

    The environment's name can be anything you liked without space in the name string.

3. install rust toolchain

    Then let's follow the [official guidance](https://www.rust-lang.org/tools/install) to install the latest Rust toolchain. This should be done by typing the following commands.

    ```bash
     $> curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
     $> rustup update
    ```

4. compile and install crystallus

    This is the last step.
    We use [pyo3](https://github.com/PyO3/pyo3) to provide you the python bindings and use [maturin](https://github.com/PyO3/maturin) to build our rust codes into the native dynamic library.

    > Please confirm that you have installed the newest version of maturin.

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

    if you just want to test the package quickly, please use [maturin develop](https://maturin.rs/tutorial.html#build-and-install-the-module-with-maturin-develop) command.

    ```bash
    $> maturin develop
    ```

## Usage

Please see examples.
