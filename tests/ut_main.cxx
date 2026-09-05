/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_HAS_MPI
#include <mpi.h>
#endif
#ifdef GAUXC_HAS_CUDA
#include <cuda_runtime.h>
#endif

int main( int argc, char* argv[] ) {
#ifdef GAUXC_HAS_MPI
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef GAUXC_HAS_CUDA
  cudaSetDevice(rank);
#elif defined(GAUXC_HAS_SYCL)
  // SYCL has no per-process cudaSetDevice equivalent: device visibility for
  // rank-to-GPU binding is controlled externally via ZE_AFFINITY_MASK /
  // ONEAPI_DEVICE_SELECTOR, so there is nothing to do here.
#endif
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
#else
  int result = Catch::Session().run( argc, argv );
#endif
  return result;
}
