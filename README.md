# About

GauXC Copyright (c) 2020, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.

# Synopsis

GauXC is a modern, modular C++ library for the evaluation of quantities related
to the exchange-correlation (XC) energy (e.g. potential, etc) in the Gaussian
basis set discretization of Kohn-Sham density function theory (KS-DFT). GauXC
provides efficient, scalable distributed memory XC integrators for both CPU and
accelerator-based (GPU) architectures. Currently, GPU support is provided through
the [CUDA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html)
framework and is only accessible to NVIDIA GPUs. However, work is underway to
provide portable interfaces to these implementations as to allow for use with
emerging exascale and post-exascale compute resources (See e.g. PRs for
[SYCL/DPC++](https://github.com/wavefunction91/GauXC/pull/4) and 
[HIP](https://github.com/wavefunction91/GauXC/pull/5) ports of GauXC). Evaluation
of the XC functional CPU/accelerator architectures is provided by the
[ExchCXX](https://github.com/wavefunction91/ExchCXX) library. Quadratures generated
by the [IntegratorXX](https://github.com/wavefunction91/IntegratorXX).

GauXC is a work in progress. Its development has been funded by the U.S.
Department of Energy Exascale Computing Project 
([NWChemEx](https://github.com/NWChemEx-Project)).

# Design Goals

* Provide a stable, portable and high-performance implementation of numerical
integrators optimized for the evaluation of XC related quantities in Gaussian
basis set KS-DFT on CPU and accelerator based architectures.
* Develop a modern, modular, extensible C++ software infrastructure which allows
for flexible and agile development in the field of KS-DFT.

# Dependencies

* CMake (3.17+)
* BLAS / LAPACK
* [ExchCXX](https://github.com/wavefunction91/ExchCXX)
* [IntegratorXX](https://github.com/wavefunction91/IntegratorXX)
* [Gau2Grid](https://github.com/dgasmith/gau2grid)
* MPI (Optional)
* [Cereal](https://github.com/USCiLab/cereal) (Optional)
* [Eigen3](https://eigen.tuxfamily.org/dox/) (Testing Only)
* [CUDA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html)/[cuBLAS](https://docs.nvidia.com/cuda/cublas/index.html) (Required only if CUDA enabled)
* [MAGMA](https://icl.utk.edu/magma/) (Optional if CUDA enabled)



# Build Instructions

GauXC provides a CMake build system with automatic dependency management (through [FetchContent]()).
As such, a simple CMake invocation will often suffice for most purposes
```
cmake -H/path/to/gauxc -B/path/to/build [GauXC configure options]
cmake --build /path/to/build
```

Currently, GauXC is not installable as a separate library (Coming Soon). GauXC may be integrated 
into existing CMake projects via `FetchContent` with the following
```
include(FetchContent)

set( GAUXC_ENABLE_TESTS OFF CACHE BOOL "Disable GauXC Tests" )
# Other GauXC Options (see below)

FetchContent_Declare( gauxc GIT_REPOSITORY https://github/com/wavefunction91/GauXC.git GIT_TAG master )
FetchContent_MakeAvailable( gauxc )


# Create my_target
target_link_libraries( my_target PUBLIC gauxc )
```

## Influential CMake Variables

| Variable Name         | Description                                                    | Default  |
|-----------------------|----------------------------------------------------------------|----------|
| `GAUXC_ENABLE_TESTS`  | Enable Testing Framework (Catch2)                              | `ON`     |
| `GAUXC_ENABLE_CEREAL` | Enable Cereal serialization  (required for testing)            | `ON`     |
| `GAUXC_ENABLE_CUDA`   | Enable CUDA XC integrator                                      | `OFF`    |
| `GAUXC_ENABLE_MAGMA`  | Enable MAGMA for batched BLAS (No effect if CUDA disabled)     | `ON`     | 

# Example Usage

Coming Soon....


# License

GauXC is made freely available under the terms of a modified 3-Clause BSD license. See
LICENSE.txt for details.

# Acknowledgments

The development of GauXC is supported by the Exascale Computing Project
(17-SC-20-SC), a collaborative effort of the U.S. Department of Energy Office
of Science and the National Nuclear Security Administration.
