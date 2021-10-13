#include "ut_common.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include <Eigen/Core>

using namespace GauXC;

void test_xc_integrator( ExecutionSpace ex, GAUXC_MPI_CODE( MPI_Comm comm, ) 
  std::string reference_file, 
  size_t quad_pad_value,
  std::string integrator_kernel = "Default",  
  std::string reduction_kernel  = "Default",
  std::string lwd_kernel        = "Default" ) {

  // Read the reference file
  using matrix_type = Eigen::MatrixXd;
  Molecule mol;
  BasisSet<double> basis;
  matrix_type P, VXC_ref;
  double EXC_ref;
  {
    read_hdf5_record( mol,   reference_file, "/MOLECULE" );
    read_hdf5_record( basis, reference_file, "/BASIS"    );

    HighFive::File file( reference_file, HighFive::File::ReadOnly );
    auto dset = file.getDataSet("/DENSITY");
    auto dims = dset.getDimensions();
    P       = matrix_type( dims[0], dims[1] );
    VXC_ref = matrix_type( dims[0], dims[1] );

    dset.read( P.data() );
    dset = file.getDataSet("/VXC");
    dset.read( VXC_ref.data() );

    dset = file.getDataSet("/EXC");
    dset.read( &EXC_ref );
    
  }


  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
  //LoadBalancerFactory lb_factory(ExecutionSpace::Host, "REPLICATED-FILLIN");
  auto lb = lb_factory.get_instance(GAUXC_MPI_CODE(comm,) mol, mg, basis, 
    quad_pad_value);

  functional_type func( ExchCXX::Backend::builtin, ExchCXX::Functional::PBE0, 
    ExchCXX::Spin::Unpolarized );

  XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", 
    integrator_kernel, lwd_kernel, reduction_kernel );
  auto integrator = integrator_factory.get_instance( func, lb );


  auto [ EXC, VXC ] = integrator.eval_exc_vxc( P );
  CHECK( EXC == Approx( EXC_ref ) );

  //std::cout << "VXC" << std::endl;
  //std::cout << VXC << std::endl;
  //std::cout << "VXC_ref" << std::endl;
  //std::cout << VXC_ref << std::endl;
  //std::cout << "DIFF" << std::endl;
  //std::cout << (VXC-VXC_ref) << std::endl;

  auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
  CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 ); 

  // Check if the integrator propagates state correctly
  if( true ) {
    auto [ EXC1, VXC1 ] = integrator.eval_exc_vxc( P );
    CHECK( EXC1 == Approx( EXC_ref ) );
    auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
    CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 ); 
  }


#if 0
  if( ex == ExecutionSpace::Host ) {
    std::cout << "EXC = " << EXC << std::endl;
    auto EXC_GRAD = integrator.eval_exc_grad( P );
    for( auto x : EXC_GRAD ) std::cout << x << std::endl;
  }
#endif
}


TEST_CASE( "Benzene / PBE0 / cc-pVDZ", "[xc-integrator]" ) {

  const std::string reference_file = 
    GAUXC_REF_DATA_PATH "/benzene_pbe0_cc-pvdz_ufg_ssf.hdf5";

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host" ) {
    test_xc_integrator( ExecutionSpace::Host, GAUXC_MPI_CODE(comm,) reference_file,
      1 );
  }
#endif

#ifdef GAUXC_ENABLE_DEVICE
  SECTION( "Device" ) {
    SECTION( "Incore - MPI Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, 32, "Default" );
    }

    #ifdef GAUXC_ENABLE_MAGMA
    SECTION( "Incore - MPI Reduction - MAGMA" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, 32, "Default", "Default", "Scheme1-MAGMA" );
    }
    #endif

    #ifdef GAUXC_ENABLE_NCCL
    SECTION( "Incore - NCCL Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,)
        reference_file, 32, "Default", "NCCL" );
    }
    #endif

    SECTION( "ShellBatched" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, 32, "ShellBatched" );
    }
  }
#endif

}
