/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
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
  std::shared_ptr<functional_type> func, 
  PruningScheme pruning_scheme,
  bool check_grad,
  bool check_integrate_den,
  bool check_k,
  std::string integrator_kernel = "Default",  
  std::string reduction_kernel  = "Default",
  std::string lwd_kernel        = "Default",
  std::shared_ptr<functional_type> epcfunc = nullptr) {

  // Read the reference file
  using matrix_type = Eigen::MatrixXd;
  Molecule mol;
  BasisSet<double> basis;
  matrix_type P, Pz, Py, Px, VXC_ref, VXCz_ref, VXCy_ref, VXCx_ref, K_ref;
  double EXC_ref;
  std::vector<double> EXC_GRAD_ref;
  bool has_k = false, has_exc_grad = false, rks = true, uks = false, gks = false;
  // For NEO-DFT Test:
  bool neo = false;
  BasisSet<double> protonic_basis;
  matrix_type protonic_Ps, protonic_Pz, protonic_VXCs_ref, protonic_VXCz_ref;
  double protonic_EXC_ref;
  {
    read_hdf5_record( mol,   reference_file, "/MOLECULE" );
    read_hdf5_record( basis, reference_file, "/BASIS"    );

    HighFive::File file( reference_file, HighFive::File::ReadOnly );
    
    std::string den="/DENSITY";
    std::string den2="/DENSITY_Z";
    std::string den3="/DENSITY_Y";
    std::string den4="/DENSITY_X";
    std::string vxc="/VXC";
    std::string vxc2="VXC_Z";
    std::string vxc3="VXC_Y";
    std::string vxc4="VXC_X";

    if (file.exist("/DENSITY_Z")) { rks = false; }

    if (file.exist("/DENSITY_Z") and not file.exist("/DENSITY_Y") and not file.exist("/DENSITY_X")) {
       den="/DENSITY_SCALAR";
       vxc="/VXC_SCALAR";
       uks=true;
    }
     
    if (file.exist("/DENSITY_X") and file.exist("/DENSITY_Y") and file.exist("/DENSITY_Z")) {
       den="/DENSITY_SCALAR";
       vxc="/VXC_SCALAR";
       gks=true;
    }
 
    if (file.exist("/PROTONIC_DENSITY_SCALAR") and file.exist("/PROTONIC_DENSITY_Z")) 
       neo=true;
    
    if(neo and !integrator_kernel.compare("ShellBatched")) return;
 
    auto dset = file.getDataSet(den);
    
    auto dims = dset.getDimensions();

    P        = matrix_type( dims[0], dims[1] );
    VXC_ref  = matrix_type( dims[0], dims[1] );
    if (not rks) {
      Pz       = matrix_type( dims[0], dims[1] );
      VXCz_ref = matrix_type( dims[0], dims[1] );
    } 
    if (gks) {
      Py       = matrix_type( dims[0], dims[1] );
      VXCy_ref = matrix_type( dims[0], dims[1] );
      Px       = matrix_type( dims[0], dims[1] );
      VXCx_ref = matrix_type( dims[0], dims[1] );
    }


    dset.read( P.data() );
    dset = file.getDataSet(vxc);
    dset.read( VXC_ref.data() );

    if (not rks) {
      dset = file.getDataSet(den2);
      dset.read( Pz.data() );
      dset = file.getDataSet(vxc2);
      dset.read( VXCz_ref.data() );
    }

    if (gks) {
      dset = file.getDataSet(den3);
      dset.read( Py.data() );
      dset = file.getDataSet(vxc3);
      dset.read( VXCy_ref.data() );
      dset = file.getDataSet(den4);
      dset.read( Px.data() );
      dset = file.getDataSet(vxc4);
      dset.read( VXCx_ref.data() );
    }    

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

    if (neo) {

      read_hdf5_record( protonic_basis, reference_file, "/PROTONIC_BASIS"    );
          
      std::string prot_den_s="/PROTONIC_DENSITY_SCALAR";
      std::string prot_den_z="/PROTONIC_DENSITY_Z";
      std::string prot_vxc_s="/PROTONIC_VXC_SCALAR";
      std::string prot_vxc_z="/PROTONIC_VXC_Z";

      auto protonic_dset = file.getDataSet(prot_den_s);
      auto protonic_dims = protonic_dset.getDimensions();

      protonic_Ps        = matrix_type( protonic_dims[0], protonic_dims[1] );
      protonic_Pz        = matrix_type( protonic_dims[0], protonic_dims[1] );
      protonic_VXCs_ref  = matrix_type( protonic_dims[0], protonic_dims[1] );
      protonic_VXCz_ref  = matrix_type( protonic_dims[0], protonic_dims[1] );

      protonic_dset.read( protonic_Ps.data() );
      protonic_dset = file.getDataSet(prot_vxc_s);
      protonic_dset.read( protonic_VXCs_ref.data() );

      protonic_dset = file.getDataSet(prot_den_z);
      protonic_dset.read( protonic_Pz.data() );
      protonic_dset = file.getDataSet(prot_vxc_z);
      protonic_dset.read( protonic_VXCz_ref.data() );

      protonic_dset = file.getDataSet("/PROTONIC_EXC");
      protonic_dset.read( &protonic_EXC_ref );
    }
  }

  if( (uks or gks) and ex == ExecutionSpace::Device and func->is_mgga() ) return;

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  auto mg = MolGridFactory::create_default_molgrid(mol, pruning_scheme,
    BatchSize(512), RadialQuad::MuraKnowles, AtomicGridSizeDefault::UltraFineGrid);

  // Construct Load Balancer
  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
  std::unique_ptr<LoadBalancer> lb;
  if (neo) {
    for( auto& sh : protonic_basis ) 
      sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
    lb = std::make_unique<LoadBalancer>( lb_factory.get_instance(rt, mol, mg, basis, protonic_basis) );
  } else{
    lb = std::make_unique<LoadBalancer>( lb_factory.get_instance(rt, mol, mg, basis) );
  }

  // Construct Weights Module
  MolecularWeightsFactory mw_factory( ex, "Default", MolecularWeightsSettings{} );
  auto mw = mw_factory.get_instance();

  // Apply partition weights
  mw.modify_weights(*lb);

  // Construct XC Functional
  //auto Spin = uks ? ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
  //functional_type func( ExchCXX::Backend::builtin, func_key, Spin );

  // Construct XCIntegrator
  XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", 
    integrator_kernel, lwd_kernel, reduction_kernel );
  std::unique_ptr<XCIntegrator<matrix_type>> integrator;
  if(neo){
    integrator = std::make_unique<XCIntegrator<matrix_type>>( integrator_factory.get_instance(*func, *epcfunc, *lb) );
  }else{
    integrator = std::make_unique<XCIntegrator<matrix_type>>( integrator_factory.get_instance(*func, *lb) );
  }

  // Integrate Density
  if( check_integrate_den and rks) {
    auto N_EL_ref = std::accumulate( mol.begin(), mol.end(), 0ul,
      [](const auto& a, const auto &b) { return a + b.Z.get(); });
    auto N_EL = integrator->integrate_den( P );
    // Factor of 2 b/c P is the alpha density for RKS
    CHECK( N_EL == Approx(N_EL_ref/2.0).epsilon(1e-6) );
  }

  // Integrate EXC/VXC
  if ( rks ) {
    double EXC, protonic_EXC;
    matrix_type VXC, protonic_VXCs, protonic_VXCz;

    if(neo) std::tie( EXC, protonic_EXC, VXC, protonic_VXCs, protonic_VXCz ) = integrator->neo_eval_exc_vxc( P, protonic_Ps, protonic_Pz );
    else    std::tie( EXC, VXC ) = integrator->eval_exc_vxc( P );

    // Check EXC/VXC
    auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
    CHECK( EXC == Approx( EXC_ref ) );
    CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 ); 

    if ( neo ) {
      auto protonic_VXCs_diff_nrm = ( protonic_VXCs - protonic_VXCs_ref ).norm();
      auto protonic_VXCz_diff_nrm = ( protonic_VXCz - protonic_VXCz_ref ).norm();
      CHECK( protonic_EXC == Approx( protonic_EXC_ref ) );
      CHECK( protonic_VXCs_diff_nrm / protonic_basis.nbf() < 1e-10 );
      CHECK( protonic_VXCz_diff_nrm / protonic_basis.nbf() < 1e-10 );
    }

    // Check if the integrator propagates state correctly
    { 
      double EXC1, protonic_EXC1;
      matrix_type VXC1, protonic_VXCs1, protonic_VXCz1;

      if(neo) std::tie( EXC1, protonic_EXC1, VXC1, protonic_VXCs1, protonic_VXCz1 ) = integrator->neo_eval_exc_vxc( P, protonic_Ps, protonic_Pz );
      else    std::tie( EXC1, VXC1 ) = integrator->eval_exc_vxc( P );
      
      CHECK( EXC1 == Approx( EXC_ref ) );
      auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
      CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 ); 

      if ( neo ) {
        auto protonic_VXCs1_diff_nrm = ( protonic_VXCs1 - protonic_VXCs_ref ).norm();
        auto protonic_VXCz1_diff_nrm = ( protonic_VXCz1 - protonic_VXCz_ref ).norm();
        CHECK( protonic_EXC1 == Approx( protonic_EXC_ref ) );
        CHECK( protonic_VXCs1_diff_nrm / protonic_basis.nbf() < 1e-10 );
        CHECK( protonic_VXCz1_diff_nrm / protonic_basis.nbf() < 1e-10 );
      }
    }

    // Check EXC-only path
    if(neo) return; // NEO EXC-only NYI
    auto EXC2 = integrator->eval_exc( P );
    CHECK(EXC2 == Approx(EXC));

  } else if (uks) {
    double EXC, protonic_EXC;
    matrix_type VXC, VXCz, protonic_VXCs, protonic_VXCz;

    if(neo) std::tie( EXC, protonic_EXC, VXC, VXCz, protonic_VXCs, protonic_VXCz ) = integrator->neo_eval_exc_vxc( P, Pz, protonic_Ps, protonic_Pz );
    else    std::tie( EXC, VXC, VXCz ) = integrator->eval_exc_vxc( P, Pz );

    // Check EXC/VXC
    auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
    auto VXCz_diff_nrm = ( VXCz - VXCz_ref ).norm();
    CHECK( EXC == Approx( EXC_ref ) );
    CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 );
    CHECK( VXCz_diff_nrm / basis.nbf() < 1e-10 );
    
    if ( neo ) {
      auto protonic_VXCs_diff_nrm = ( protonic_VXCs - protonic_VXCs_ref ).norm();
      auto protonic_VXCz_diff_nrm = ( protonic_VXCz - protonic_VXCz_ref ).norm();
      CHECK( protonic_EXC == Approx( protonic_EXC_ref ) );
      CHECK( protonic_VXCs_diff_nrm / protonic_basis.nbf() < 1e-10 );
      CHECK( protonic_VXCz_diff_nrm / protonic_basis.nbf() < 1e-10 );
    }

    // Check if the integrator propagates state correctly
    { 
      double EXC1, protonic_EXC1;
      matrix_type VXC1, VXCz1, protonic_VXCs1, protonic_VXCz1;

      if(neo) std::tie( EXC1, protonic_EXC1, VXC1, VXCz1, protonic_VXCs1, protonic_VXCz1 ) = integrator->neo_eval_exc_vxc( P, Pz, protonic_Ps, protonic_Pz );
      else    std::tie( EXC1, VXC1, VXCz1 ) = integrator->eval_exc_vxc( P, Pz );
      
      auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
      auto VXCz1_diff_nrm = ( VXCz1 - VXCz_ref ).norm();
      CHECK( EXC1 == Approx( EXC_ref ) );
      CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 );
      CHECK( VXCz1_diff_nrm / basis.nbf() < 1e-10 );

      if ( neo ) {
        auto protonic_VXCs1_diff_nrm = ( protonic_VXCs1 - protonic_VXCs_ref ).norm();
        auto protonic_VXCz1_diff_nrm = ( protonic_VXCz1 - protonic_VXCz_ref ).norm();
        CHECK( protonic_EXC1 == Approx( protonic_EXC_ref ) );
        CHECK( protonic_VXCs1_diff_nrm / protonic_basis.nbf() < 1e-10 );
        CHECK( protonic_VXCz1_diff_nrm / protonic_basis.nbf() < 1e-10 );
      }
    }

    // Check EXC-only path
    if(neo) return; // NEO EXC-only NYI
    auto EXC2 = integrator->eval_exc( P, Pz );
    CHECK(EXC2 == Approx(EXC));
  } else if (gks) {
    auto [ EXC, VXC, VXCz, VXCy, VXCx ] = integrator->eval_exc_vxc( P, Pz, Py, Px );

    // Check EXC/VXC
    auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
    auto VXCz_diff_nrm = ( VXCz - VXCz_ref ).norm();
    auto VXCy_diff_nrm = ( VXCy - VXCy_ref ).norm();
    auto VXCx_diff_nrm = ( VXCx - VXCx_ref ).norm();

    CHECK( EXC == Approx( EXC_ref ) );
    CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 );
    CHECK( VXCz_diff_nrm / basis.nbf() < 1e-10 );
    CHECK( VXCy_diff_nrm / basis.nbf() < 1e-10 );
    CHECK( VXCx_diff_nrm / basis.nbf() < 1e-10 );
    // Check if the integrator propagates state correctly
    {
      auto [ EXC1, VXC1, VXCz1, VXCy1, VXCx1] = integrator->eval_exc_vxc( P, Pz, Py, Px );
      CHECK( EXC1 == Approx( EXC_ref ) );
      auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
      auto VXCz1_diff_nrm = ( VXCz1 - VXCz_ref ).norm();
      auto VXCy1_diff_nrm = ( VXCy1 - VXCy_ref ).norm();
      auto VXCx1_diff_nrm = ( VXCx1 - VXCx_ref ).norm();
      CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 );
      CHECK( VXCz1_diff_nrm / basis.nbf() < 1e-10 );
      CHECK( VXCy1_diff_nrm / basis.nbf() < 1e-10 );
      CHECK( VXCx1_diff_nrm / basis.nbf() < 1e-10 );
    }

    // Check EXC-only path
    auto EXC2 = integrator->eval_exc( P, Pz, Py, Px );
    CHECK(EXC2 == Approx(EXC));
  }



  // Check EXC Grad
  if( check_grad and has_exc_grad and rks) {
    auto EXC_GRAD = integrator->eval_exc_grad( P );
    using map_type = Eigen::Map<Eigen::MatrixXd>;
    map_type EXC_GRAD_ref_map( EXC_GRAD_ref.data(), mol.size(), 3 );
    map_type EXC_GRAD_map( EXC_GRAD.data(), mol.size(), 3 );
    auto EXC_GRAD_diff_nrm = (EXC_GRAD_ref_map - EXC_GRAD_map).norm();
    CHECK( EXC_GRAD_diff_nrm / std::sqrt(3.0*mol.size()) < 1e-10 );
  }

  // Check K
  if( has_k and check_k and rks ) {
    auto max_l = basis.max_l();
    if(max_l > 2 and ex == ExecutionSpace::Device) {
      std::cout << "Skiping device sn-K + L > 2" << std::endl;
      return;
    }
    auto K = integrator->eval_exx( P );
    CHECK((K - K.transpose()).norm() < std::numeric_limits<double>::epsilon()); // Symmetric
    CHECK( (K - K_ref).norm() / basis.nbf() < 1e-7 );
  }

}

void test_integrator(std::string reference_file, std::shared_ptr<functional_type> func, PruningScheme pruning_scheme,
  std::shared_ptr<functional_type> epcfunc = nullptr) {

#ifdef GAUXC_HAS_DEVICE
  auto rt = DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
#else
  auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#endif

#ifdef GAUXC_HAS_HOST
    SECTION( "Host" ) {
      SECTION("Reference") {
        test_xc_integrator( ExecutionSpace::Host, rt, reference_file, func,
          pruning_scheme, true, true, true, "Default", "Default", "Default", epcfunc);
      }
      SECTION("ShellBatched") {
        test_xc_integrator( ExecutionSpace::Host, rt, reference_file, func,
          pruning_scheme, false, false, false, "ShellBatched", "Default", "Default", epcfunc);
      }
    }
#endif

#ifdef GAUXC_HAS_DEVICE
  SECTION( "Device" ) {
    bool check_grad = true;
    bool check_k    = true;
    #ifdef GAUXC_HAS_HIP
    check_grad = false;
    check_k    = false;
    #endif
    SECTION( "Incore - MPI Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, rt,
        reference_file, func, pruning_scheme,  
        check_grad, true, check_k, "Default" );
    }

    #ifdef GAUXC_HAS_MAGMA
    SECTION( "Incore - MPI Reduction - MAGMA" ) {
      if(not func->is_mgga() and not func->is_polarized()) {
        test_xc_integrator( ExecutionSpace::Device, rt,
          reference_file, func, pruning_scheme,
          false, true, check_k, "Default", "Default", 
          "Scheme1-MAGMA" );
      }
    }
    #endif

    #ifdef GAUXC_HAS_CUTLASS
    SECTION( "Incore - MPI Reduction - CUTLASS" ) {
      if(not func->is_mgga() and not func->is_polarized()) {
        test_xc_integrator( ExecutionSpace::Device, rt, 
          reference_file, func, pruning_scheme,
          false, true, false, "Default", "Default", 
          "Scheme1-CUTLASS" );
      }
    }
    #endif


    #ifdef GAUXC_HAS_NCCL
    SECTION( "Incore - NCCL Reduction" ) {
      test_xc_integrator( ExecutionSpace::Device, rt,
        reference_file, func, pruning_scheme, 
        false, false, false, "Default", "NCCL" );
    }
    #endif

    SECTION( "ShellBatched" ) {
      test_xc_integrator( ExecutionSpace::Device, rt, 
        reference_file, func, pruning_scheme,  
        false, false, false, "ShellBatched" );
    }
  }
#endif

}

std::shared_ptr<functional_type> make_functional(ExchCXX::Functional func_key, ExchCXX::Spin spin) {
  return std::make_shared<functional_type>(ExchCXX::Backend::builtin, func_key, spin);
}


TEST_CASE( "XC Integrator", "[xc-integrator]" ) {

  auto pol     = ExchCXX::Spin::Polarized;
  auto unpol   = ExchCXX::Spin::Unpolarized;
  auto svwn5   = ExchCXX::Functional::SVWN5;
  auto pbe0    = ExchCXX::Functional::PBE0;
  auto blyp    = ExchCXX::Functional::BLYP;
  auto scan    = ExchCXX::Functional::SCAN;
  auto r2scanl = ExchCXX::Functional::R2SCANL;
  auto epc17_2 = ExchCXX::Functional::EPC17_2;
  auto epc18_2 = ExchCXX::Functional::EPC18_2;

  // LDA Test
  SECTION( "Benzene / SVWN5 / cc-pVDZ" ) {
    auto func = make_functional(svwn5, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf.hdf5", 
        func, PruningScheme::Unpruned );
  }
  SECTION( "Benzene / SVWN5 / cc-pVDZ (Treutler)" ) {
    auto func = make_functional(svwn5, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf_treutler_prune.hdf5", 
        func, PruningScheme::Treutler );
  }
  SECTION( "Benzene / SVWN5 / cc-pVDZ (Robust)" ) {
    auto func = make_functional(svwn5, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_svwn5_cc-pvdz_ufg_ssf_robust_prune.hdf5", 
        func, PruningScheme::Robust );
  }

  // GGA Test
  SECTION( "Benzene / PBE0 / cc-pVDZ" ) {
    auto func = make_functional(pbe0, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_pbe0_cc-pvdz_ufg_ssf.hdf5", 
        func, PruningScheme::Unpruned );
  }

  // MGGA Test (TAU Only)
  SECTION( "Cytosine / SCAN / cc-pVDZ") {
    auto func = make_functional(scan, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/cytosine_scan_cc-pvdz_ufg_ssf_robust.hdf5", 
        func, PruningScheme::Robust );
  }

  // MGGA Test (TAU + LAPL)
  SECTION( "Cytosine / R2SCANL / cc-pVDZ") {
    auto func = make_functional(r2scanl, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/cytosine_r2scanl_cc-pvdz_ufg_ssf_robust.hdf5", 
        func, PruningScheme::Robust );
  }

  //UKS LDA Test
  SECTION( "Li / SVWN5 / sto-3g" ) {
    auto func = make_functional(svwn5, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/li_svwn5_sto3g_uks.bin",
        func, PruningScheme::Unpruned );
  }

  //UKS GGA Test
  SECTION( "Li / BLYP / sto-3g" ) {
    auto func = make_functional(blyp, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/li_blyp_sto3g_uks.bin",
        func, PruningScheme::Unpruned );
  }

  // UKS MGGA Test (TAU Only)
  SECTION( "Cytosine (doublet) / SCAN / cc-pVDZ") {
    auto func = make_functional(scan, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/cytosine_scan_cc-pvdz_ufg_ssf_robust_uks.hdf5", 
        func, PruningScheme::Robust );
  }

  // UKS MGGA Test (TAU + LAPL)
  SECTION( "Cytosine (doublet) / R2SCANL / cc-pVDZ") {
    auto func = make_functional(r2scanl, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/cytosine_r2scanl_cc-pvdz_ufg_ssf_robust_uks.hdf5", 
        func, PruningScheme::Robust );
  }

  // GKS GGA Test
  SECTION( "H3 / BLYP / cc-pvdz" ) {
    auto func = make_functional(blyp, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/h3_blyp_cc-pvdz_ssf_gks.bin",
        func, PruningScheme::Unpruned );
  }

  // sn-LinK Test
  SECTION( "Benzene / PBE0 / 6-31G(d)" ) {
    auto func = make_functional(pbe0, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/benzene_631gd_pbe0_ufg.hdf5", 
        func, PruningScheme::Unpruned );
  }

  // sn-LinK + f functions
  SECTION( "H2O2 / PBE0 / def2-TZVP" ) {
    auto func = make_functional(pbe0, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/h2o2_def2-tzvp.hdf5", 
        func, PruningScheme::Unpruned );
  }

  // sn-LinK + g functions
  SECTION( "H2O2 / PBE0 / def2-QZVP" ) {
    auto func = make_functional(pbe0, unpol);
    test_integrator(GAUXC_REF_DATA_PATH "/h2o2_def2-qzvp.hdf5", 
        func, PruningScheme::Unpruned );
  }

  // EPC Tests
  // epc-17-2 Test (small basis)
  SECTION( "COH2 / BLYP,EPC-17-2 / sto-3g, prot-sp" ) {
    auto func = make_functional(blyp, unpol);
    auto epcfunc = make_functional(epc17_2, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/coh2_blyp_epc17-2_sto-3g_protsp_ssf.hdf5", 
        func, PruningScheme::Unpruned, epcfunc);
  }
  // epc-17-2 Test (larger basis)
  SECTION( "COH2 / BLYP,EPC-17-2 / cc-pVDZ, prot-PB4-D" ) {
    auto func = make_functional(blyp, unpol);
    auto epcfunc = make_functional(epc17_2, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/coh2_blyp_epc17-2_cc-pvdz_pb4d_ssf.hdf5", 
        func, PruningScheme::Unpruned, epcfunc);
  }
  // epc-18-2 Test
  SECTION( "COH2 / BLYP,EPC-18-2 / cc-pVDZ, prot-PB4-D" ) {
    auto func = make_functional(blyp, unpol);
    auto epcfunc = make_functional(epc18_2, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/coh2_blyp_epc18-2_cc-pvdz_pb4d_ssf.hdf5", 
        func, PruningScheme::Unpruned, epcfunc);
   }
  // UKS NEO epc-17-2 Test 
  SECTION( "OH2+ / BLYP,EPC-17-2 / cc-pVDZ, prot-PB4-D" ) {
    auto func = make_functional(blyp, pol);
    auto epcfunc = make_functional(epc17_2, pol);
    test_integrator(GAUXC_REF_DATA_PATH "/coh2_blyp_epc17-2_cc-pvdz_pb4d_ssf_uks.hdf5", 
        func, PruningScheme::Unpruned, epcfunc);
  }
} 
  
