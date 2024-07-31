/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_HAS_MPI
#include <mpi.h>
#endif

int main( int argc, char* argv[] ) {
#ifdef GAUXC_HAS_MPI
  MPI_Init(&argc, &argv);
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
#else
  int result = Catch::Session().run( argc, argv );
#endif
  return result;
}
