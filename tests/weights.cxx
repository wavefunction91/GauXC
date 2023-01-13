#include "ut_common.hpp"
#include <gauxc/molgrid.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/load_balancer.hpp>
#include <gauxc/util/div_ceil.hpp>
#include <fstream>
#include <string>

#include "weights_generate.hpp"
#include "weights_host.hpp"
#include "weights_cuda.hpp"
#include "weights_hip.hpp"

//#define GENERATE_TESTS
TEST_CASE( "Benzene", "[weights]" ) {

  auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#ifdef GENERATE_TESTS
  if(rt.comm_size() > 1) return;
#endif

  Molecule mol = make_benzene();

#ifdef GENERATE_TESTS
  BasisSet<double> basis = make_631Gd( mol, SphericalType(true) );
  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

  std::ofstream ref_data( "benzene_weights_ssf.bin", std::ios::binary );
  generate_weights_data( mol, basis, ref_data );  
#else

  std::ifstream ref_data( GAUXC_REF_DATA_PATH "/benzene_weights_ssf.bin", 
                          std::ios::binary );

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host Weights" ) {
    test_host_weights( ref_data );
  }
#endif

#ifdef GAUXC_ENABLE_DEVICE
  SECTION( "Device Weights" ) {
#ifdef GAUXC_ENABLE_CUDA
    test_cuda_weights( ref_data );
#elif defined(GAUXC_ENABLE_HIP)
    test_hip_weights( ref_data );
#endif
  }
#endif
#endif
}


