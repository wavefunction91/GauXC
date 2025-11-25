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
#include "collocation_common.hpp"
#include "collocation_host.hpp"
#include "collocation_cuda.hpp"
#include "collocation_hip.hpp"

//#define GENERATE_TESTS

#if defined(GENERATE_TESTS) && !defined(GAUXC_HAS_HOST)
  #error "Host Integrator Must Be Enabled to Generate Tests"
#endif

TEST_CASE( "Water / cc-pVDZ", "[collocation]" ) {

#ifdef GENERATE_TESTS
#ifdef GAUXC_HAS_MPI
  int world_size;
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  if( world_size > 1 ) return;
#endif
#endif

  Molecule mol           = make_water();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

#ifdef GENERATE_TESTS

  generate_collocation_data( mol, basis, "water_cc-pVDZ_collocation.hdf5" );

#else

  std::string ref_file = GAUXC_REF_DATA_PATH "/water_cc-pVDZ_collocation.hdf5";

#ifdef GAUXC_HAS_HOST
  SECTION( "Host Eval" ) {
    test_host_collocation( basis, ref_file );
  }

  SECTION( "Host Eval Grad" ) {
    test_host_collocation_deriv1( basis, ref_file );
  }

  SECTION( "Host Eval Hessian" ) {
    test_host_collocation_deriv2( basis, ref_file );
  }
#endif

#ifdef GAUXC_HAS_CUDA
  BasisSetMap basis_map( basis, mol );
  SECTION( "CUDA Eval" ) {
    test_cuda_collocation_masked_combined( basis, ref_file, false );
  }

  SECTION( "CUDA Eval Grad" ) {
    test_cuda_collocation_masked_combined( basis, ref_file, true );
  }

  SECTION( "CUDA Eval Hessian" ) {
    test_cuda_collocation_masked_combined_deriv2( basis, ref_file, false, false );
  }

  SECTION( "CUDA Eval Laplacian" ) {
    test_cuda_collocation_masked_combined_deriv2( basis, ref_file, true, false );
  }

  SECTION( "CUDA Eval Laplacian Gradient" ) {
    test_cuda_collocation_masked_combined_deriv2( basis, ref_file, true, true );
  }
#endif // GAUXC_HAS_CUDA

#ifdef GAUXC_HAS_HIP
  SECTION( "HIP Eval" ) {
    test_hip_collocation_masked_combined( basis, ref_file, false );
  }

  SECTION( "HIP Eval Grad" ) {
    test_hip_collocation_masked_combined( basis, ref_file, true );
  }
#endif // GAUXC_HAS_HIP




#endif

}
