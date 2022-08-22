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
  ExchCXX::Functional func_key, 
  size_t quad_pad_value,
  bool check_grad,
  std::string integrator_kernel = "Default",  
  std::string reduction_kernel  = "Default",
  std::string lwd_kernel        = "Default" ) {

  // Read the reference file
  using matrix_type = Eigen::MatrixXd;
  Molecule mol;
  BasisSet<double> basis;
  matrix_type P, VXC_ref;
  double EXC_ref;
  std::vector<double> EXC_GRAD_ref;
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

    EXC_GRAD_ref.resize( 3*mol.size() );
    dset = file.getDataSet("/EXC_GRAD");
    dset.read( EXC_GRAD_ref.data() );
    
  }


  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
  //LoadBalancerFactory lb_factory(ExecutionSpace::Host, "REPLICATED-FILLIN");
  auto lb = lb_factory.get_instance(GAUXC_MPI_CODE(comm,) mol, mg, basis, 
    quad_pad_value);

  functional_type func( ExchCXX::Backend::builtin, func_key, ExchCXX::Spin::Unpolarized );

  XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", 
    integrator_kernel, lwd_kernel, reduction_kernel );
  auto integrator = integrator_factory.get_instance( func, lb );

  //auto N_EL = integrator.integrate_den( P );

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


  // Check EXC Grad
  if( check_grad ) {
    auto EXC_GRAD = integrator.eval_exc_grad( P );
    using map_type = Eigen::Map<Eigen::MatrixXd>;
    map_type EXC_GRAD_ref_map( EXC_GRAD_ref.data(), mol.size(), 3 );
    map_type EXC_GRAD_map( EXC_GRAD.data(), mol.size(), 3 );
    auto EXC_GRAD_diff_nrm = (EXC_GRAD_ref_map - EXC_GRAD_map).norm();
    CHECK( EXC_GRAD_diff_nrm / std::sqrt(3.0*mol.size()) < 1e-10 );
  }

}

void test_integrator(std::string reference_file, ExchCXX::Functional func) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host" ) {
    test_xc_integrator( ExecutionSpace::Host, GAUXC_MPI_CODE(comm,) reference_file, func,
      1, true );
  }
#endif

#ifdef GAUXC_ENABLE_DEVICE
  SECTION( "Device" ) {
    bool check_grad = true;
    #ifdef GAUXC_ENABLE_HIP
    check_grad = false;
    #endif
    SECTION( "Incore - MPI Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, func, 32, check_grad, "Default" );
    }

    #ifdef GAUXC_ENABLE_MAGMA
    SECTION( "Incore - MPI Reduction - MAGMA" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, func, 32, false, "Default", "Default", "Scheme1-MAGMA" );
    }
    #endif

    #ifdef GAUXC_ENABLE_NCCL
    SECTION( "Incore - NCCL Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,)
        reference_file, func, 32, false, "Default", "NCCL" );
    }
    #endif

    SECTION( "ShellBatched" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_CODE(comm,) 
        reference_file, func, 32, false, "ShellBatched" );
    }
  }
#endif

}


TEST_CASE( "XC Integrator", "[xc-integrator]" ) {


  // LDA Test
  SECTION( "Benzene / SVWN5 / cc-pVDZ" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf.hdf5", ExchCXX::Functional::SVWN5 );
  }

  // GGA Test
  SECTION( "Benzene / PBE0 / cc-pVDZ" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_pbe0_cc-pvdz_ufg_ssf.hdf5", ExchCXX::Functional::PBE0 );
  }

}
