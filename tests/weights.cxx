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
TEST_CASE( "Partition Weights", "[weights]" ) {

  auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#ifdef GENERATE_TESTS
  if(rt.comm_size() > 1) return;
#endif

  Molecule mol = make_benzene();

#ifdef GENERATE_TESTS
  BasisSet<double> basis = make_631Gd( mol, SphericalType(true) );
  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

  generate_weights_data( mol, basis, "benzene_weights_becke.hdf5", XCWeightAlg::Becke );  
  generate_weights_data( mol, basis, "benzene_weights_ssf.hdf5", XCWeightAlg::SSF );  
  generate_weights_data( mol, basis, "benzene_weights_lko.hdf5", XCWeightAlg::LKO );
  return;
#else


#ifdef GAUXC_HAS_HOST
  SECTION("Becke") {
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_becke.hdf5";
  test_host_weights( ref_file, XCWeightAlg::Becke );
  }
  SECTION("LKO") {
  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_lko.hdf5";
  test_host_weights( ref_file, XCWeightAlg::LKO );
  }
#endif


  SECTION("SSF") {

  std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_weights_ssf.hdf5";

#ifdef GAUXC_HAS_HOST
  SECTION( "Host Weights" ) {
    test_host_weights( ref_file, XCWeightAlg::SSF );
  }
#endif

#ifdef GAUXC_HAS_DEVICE
  SECTION( "Device Weights" ) {
#ifdef GAUXC_HAS_CUDA
    test_cuda_weights( ref_file );
#elif defined(GAUXC_HAS_HIP)
    test_hip_weights( ref_file );
#endif
  }
#endif
#endif

  }
}


