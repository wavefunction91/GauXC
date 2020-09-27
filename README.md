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
[SYCL](https://github.com/wavefunction91/GauXC/pull/4) and 
[HIP](https://github.com/wavefunction91/GauXC/pull/5) ports of GauXC). Evaluation
of the XC functional CPU/accelerator architectures is provided by the
[ExchCXX](https://github.com/wavefunction91/ExchCXX) library.

GauXC is a work in progress. Its development has been funded by the U.S.
Department of Energy Exascale Computing Project 
([NWChemEx](https://github.com/NWChemEx-Project)).

# Design Goals

* Provide a stable, portable and high-performance implementation of numerical
integrators optimized for the evaluation of XC related quantities in Gaussian
basis set KS-DFT on CPU and accelerator based architectures.
* Develop a modern, modular, extensible C++ software infrastructure which allows
for flexible and agile development in the field of KS-DFT.

# Example Usage

Coming Soon....


# License

GauXC is made freely available under the terms of a modified 3-Clause BSD license. See
LICENSE.txt for details.

# Acknowledgments

The development of GauXC is supported by the Exascale Computing Project
(17-SC-20-SC), a collaborative effort of the U.S. Department of Energy Office
of Science and the National Nuclear Security Administration.
