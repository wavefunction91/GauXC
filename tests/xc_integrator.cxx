/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "ut_common.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/molecular_weights.hpp>

#include <gauxc/molgrid/defaults.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include <Eigen/Core>

using namespace GauXC;

void test_xc_integrator( ExecutionSpace ex, const RuntimeEnvironment& rt,
  std::string reference_file, 
  ExchCXX::Functional func_key, 
  PruningScheme pruning_scheme,
  size_t quad_pad_value,
  bool check_grad,
  bool check_integrate_den,
  bool check_k,
  std::string integrator_kernel = "Default",  
  std::string reduction_kernel  = "Default",
  std::string lwd_kernel        = "Default" ) {

  // Read the reference file
  using matrix_type = Eigen::MatrixXd;
  Molecule mol;
  BasisSet<double> basis;
  matrix_type P, VXC_ref, K_ref;
  double EXC_ref;
  std::vector<double> EXC_GRAD_ref;
  bool has_k = false, has_exc_grad = false;
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

    has_exc_grad = file.exist("/EXC_GRAD");
    if( has_exc_grad ) {
      EXC_GRAD_ref.resize( 3*mol.size() );
      dset = file.getDataSet("/EXC_GRAD");
      dset.read( EXC_GRAD_ref.data() );
    }
    
    has_k = file.exist("/K");
    if(has_k) {
        K_ref = matrix_type(dims[0], dims[1]);
        dset = file.getDataSet("/K");
        dset.read( K_ref.data() );
    }
  }


  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  auto mg = MolGridFactory::create_default_molgrid(mol, pruning_scheme,
    BatchSize(512), RadialQuad::MuraKnowles, AtomicGridSizeDefault::UltraFineGrid);

  // Construct Load Balancer
  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
  auto lb = lb_factory.get_instance(rt, mol, mg, basis, quad_pad_value);

  // Construct Weights Module
  MolecularWeightsFactory mw_factory( ex, "Default", MolecularWeightsSettings{} );
  auto mw = mw_factory.get_instance();

  // Apply partition weights
  mw.modify_weights(lb);

  // Construct XC Functional
  functional_type func( ExchCXX::Backend::builtin, func_key, ExchCXX::Spin::Unpolarized );

  // Construct XCIntegrator
  XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", 
    integrator_kernel, lwd_kernel, reduction_kernel );
  auto integrator = integrator_factory.get_instance( func, lb );

  // Integrate Density
  if( check_integrate_den ) {
    auto N_EL_ref = std::accumulate( mol.begin(), mol.end(), 0ul,
      [](const auto& a, const auto &b) { return a + b.Z.get(); });
    auto N_EL = integrator.integrate_den( P );
    CHECK( N_EL == Approx(N_EL_ref).epsilon(1e-6) );
  }

  // Integrate EXC/VXC
  auto [ EXC, VXC ] = integrator.eval_exc_vxc( P );

  // Check EXC/VXC
  auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
  CHECK( EXC == Approx( EXC_ref ) );
  CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 ); 

  // Check if the integrator propagates state correctly
  {
    auto [ EXC1, VXC1 ] = integrator.eval_exc_vxc( P );
    CHECK( EXC1 == Approx( EXC_ref ) );
    auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
    CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 ); 
  }


  // Check EXC Grad
  if( check_grad and has_exc_grad ) {
    auto EXC_GRAD = integrator.eval_exc_grad( P );
    using map_type = Eigen::Map<Eigen::MatrixXd>;
    map_type EXC_GRAD_ref_map( EXC_GRAD_ref.data(), mol.size(), 3 );
    map_type EXC_GRAD_map( EXC_GRAD.data(), mol.size(), 3 );
    auto EXC_GRAD_diff_nrm = (EXC_GRAD_ref_map - EXC_GRAD_map).norm();
    CHECK( EXC_GRAD_diff_nrm / std::sqrt(3.0*mol.size()) < 1e-10 );
  }

  // Check K
  if( has_k and check_k ) {
    auto K = integrator.eval_exx( P );
    CHECK((K - K.transpose()).norm() < std::numeric_limits<double>::epsilon()); // Symmetric
    CHECK( (K - K_ref).norm() / basis.nbf() < 1e-7 );
  }

}

void test_integrator(std::string reference_file, ExchCXX::Functional func, PruningScheme pruning_scheme) {

#ifdef GAUXC_ENABLE_DEVICE
  auto rt = DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
#else
  auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#endif

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host" ) {
    test_xc_integrator( ExecutionSpace::Host, rt, reference_file, func,
      pruning_scheme, 1, true, true, true );
  }
#endif

#ifdef GAUXC_ENABLE_DEVICE
  SECTION( "Device" ) {
    bool check_grad = true;
    bool check_k    = true;
    #ifdef GAUXC_ENABLE_HIP
    check_grad = false;
    check_k    = false;
    #endif
    SECTION( "Incore - MPI Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, rt,
        reference_file, func, pruning_scheme, 1, 
        check_grad, true, check_k, "Default" );
    }

    #ifdef GAUXC_ENABLE_MAGMA
    SECTION( "Incore - MPI Reduction - MAGMA" ) {
      test_xc_integrator( ExecutionSpace::Device, rt,
        reference_file, func, pruning_scheme,
        1, false, true, check_k, "Default", "Default", 
        "Scheme1-MAGMA" );
    }
    #endif

    #ifdef GAUXC_ENABLE_CUTLASS
    SECTION( "Incore - MPI Reduction - CUTLASS" ) {
      test_xc_integrator( ExecutionSpace::Device, rt, 
        reference_file, func, pruning_scheme,
        1, false, true, false, "Default", "Default", 
        "Scheme1-CUTLASS" );
    }
    #endif


    #ifdef GAUXC_ENABLE_NCCL
    SECTION( "Incore - NCCL Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, rt,
        reference_file, func, pruning_scheme, 
        1, false, false, false, "Default", "NCCL" );
    }
    #endif

    SECTION( "ShellBatched" ) {
      test_xc_integrator( ExecutionSpace::Device, rt, 
        reference_file, func, pruning_scheme, 1, 
        false, false, false, "ShellBatched" );
    }
  }
#endif

}


TEST_CASE( "XC Integrator", "[xc-integrator]" ) {


  // LDA Test
  SECTION( "Benzene / SVWN5 / cc-pVDZ" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf.hdf5", 
        ExchCXX::Functional::SVWN5, PruningScheme::Unpruned );
  }
  SECTION( "Benzene / SVWN5 / cc-pVDZ (Treutler)" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf_treutler_prune.hdf5", 
        ExchCXX::Functional::SVWN5, PruningScheme::Treutler );
  }
  SECTION( "Benzene / SVWN5 / cc-pVDZ (Robust)" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf_robust_prune.hdf5", 
        ExchCXX::Functional::SVWN5, PruningScheme::Robust );
  }

  // GGA Test
  SECTION( "Benzene / PBE0 / cc-pVDZ" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_pbe0_cc-pvdz_ufg_ssf.hdf5", 
        ExchCXX::Functional::PBE0, PruningScheme::Unpruned );
  }

  // sn-LinK Test
  SECTION( "Benzene / PBE0 / 6-31G(d)" ) {
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_631gd_pbe0_ufg.hdf5", 
        ExchCXX::Functional::PBE0, PruningScheme::Unpruned );
  }
}
