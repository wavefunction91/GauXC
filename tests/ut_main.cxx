#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_ENABLE_MPI
#include <mpi.h>
#endif

int main( int argc, char* argv[] ) {
#ifdef GAUXC_ENABLE_MPI
  MPI_Init(&argc, &argv);
  int result = Catch::Session().run( argc, argv );
  MPI_Finalize();
#else
  int result = Catch::Session().run( argc, argv );
#endif
  return result;
}
