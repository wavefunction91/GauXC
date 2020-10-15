#include "collocation_common.hpp"
#include "collocation_host.hpp"
#include "collocation_cuda.hpp"
#include "collocation_sycl.hpp"

//#define GENERATE_TESTS

#if defined(GENERATE_TESTS) && !defined(GAUXC_ENABLE_HOST)
  #error "Host Integrator Must Be Enabled to Generate Tests"
#endif

TEST_CASE( "Water / cc-pVDZ", "[collocation]" ) {

#ifdef GENERATE_TESTS
#ifdef GAUXC_ENABLE_MPI
  int world_size;
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  if( world_size > 1 ) return;
#endif
#endif

  Molecule mol           = make_water();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) sh.set_shell_tolerance( 1e-6 );

#ifdef GENERATE_TESTS

  std::ofstream ref_data( "water_cc-pVDZ_collocation.bin", std::ios::binary );
  generate_collocation_data( mol, basis, ref_data );

#else

  std::ifstream ref_data( GAUXC_REF_DATA_PATH "/water_cc-pVDZ_collocation.bin",
                          std::ios::binary );

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host Eval" ) {
    test_host_collocation( basis, ref_data );
  }

  SECTION( "Host Eval Grad" ) {
    test_host_collocation_deriv1( basis, ref_data );
  }
#endif

#ifdef GAUXC_ENABLE_CUDA
  SECTION( "CUDA Eval: Petite Shell List" ) {
    test_cuda_collocation_petite( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked" ) {
    test_cuda_collocation_masked( basis, ref_data );
  }
  SECTION( "CUDA Eval: Petite Combined" ) {
    test_cuda_collocation_petite_combined( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked Combined" ) {
    test_cuda_collocation_masked_combined( basis, ref_data );
  }

  SECTION( "CUDA Eval Grad: Petite Shell List" ) {
    test_cuda_collocation_deriv1_petite( basis, ref_data );
  }
  SECTION( "CUDA Eval Grad: Masked" ) {
    test_cuda_collocation_deriv1_masked( basis, ref_data );
  }
  SECTION( "CUDA Eval Grad: Petite Combined" ) {
    test_cuda_collocation_petite_combined_deriv1( basis, ref_data );
  }
  SECTION( "CUDA Eval: Masked Combined" ) {
    test_cuda_collocation_masked_combined_deriv1( basis, ref_data );
  }
#endif // GAUXC_ENABLE_CUDA



#ifdef GAUXC_ENABLE_SYCL
//  cl::sycl::gpu_selector device_selector;
//  cl::sycl::queue syclQueue = cl::sycl::queue(device_selector,
//                                              cl::sycl::property_list{cl::sycl::property::queue::in_order{}});
//
//  SECTION( "SYCL Eval: Petite Shell List" ) {
//    test_sycl_collocation_petite( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval: Masked" ) {
//    test_sycl_collocation_masked( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval: Petite Combined" ) {
//    test_sycl_collocation_petite_combined( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval: Masked Combined" ) {
//    test_sycl_collocation_masked_combined( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//
//  SECTION( "SYCL Eval Grad: Petite Shell List" ) {
//    test_sycl_collocation_deriv1_petite( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval Grad: Masked" ) {
//    test_sycl_collocation_deriv1_masked( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval Grad: Petite Combined" ) {
//    test_sycl_collocation_petite_combined_deriv1( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
//  SECTION( "SYCL Eval: Masked Combined" ) {
//    test_sycl_collocation_masked_combined_deriv1( basis, ref_data, syclQueue );
//  }
//  syclQueue.wait_and_throw();
#endif // GAUXC_ENABLE_SYCL

#endif

}
