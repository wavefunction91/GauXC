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
#include <gauxc/util/mpi.hpp>
#include <highfive/H5File.hpp>
#include <Eigen/Core>

using namespace GauXC;


void test_onedft_integrator( ExecutionSpace ex, const RuntimeEnvironment& rt,
    std::string reference_file, 
    std::string onedft_model_path,
    std::string integrator_kernel = "Default",  
    std::string reduction_kernel  = "Default",
    std::string lwd_kernel        = "Default") {

    using matrix_type = Eigen::MatrixXd;
    Molecule mol;
    BasisSet<double> basis;
    matrix_type P, Pz, VXC_ref, VXCz_ref;
    double EXC_ref;

    read_hdf5_record( mol,   reference_file, "/MOLECULE" );
    read_hdf5_record( basis, reference_file, "/BASIS"    );

    HighFive::File file( reference_file, HighFive::File::ReadOnly );

    std::string den="/DENSITY_SCALAR";
    std::string den2="/DENSITY_Z";
    std::string vxc="/VXC_SCALAR";
    std::string vxc2="/VXC_Z";

    auto dset = file.getDataSet(den);
    auto dims = dset.getDimensions();
    P        = matrix_type( dims[0], dims[1] );
    VXC_ref  = matrix_type( dims[0], dims[1] );
    Pz       = matrix_type( dims[0], dims[1] );
    VXCz_ref = matrix_type( dims[0], dims[1] );

    dset.read( P.data() );
    dset = file.getDataSet(vxc);
    dset.read( VXC_ref.data() );
    dset = file.getDataSet(den2);
    dset.read( Pz.data() );
    dset = file.getDataSet(vxc2);
    dset.read( VXCz_ref.data() );

    dset = file.getDataSet("/EXC");
    dset.read( &EXC_ref );

    auto mg = MolGridFactory::create_default_molgrid(mol, PruningScheme::Unpruned,
        BatchSize(512), RadialQuad::MuraKnowles, AtomicGridSizeDefault::UltraFineGrid);

    LoadBalancerFactory lb_factory(ex, "Default");
    auto lb = lb_factory.get_instance(rt, mol, mg, basis);
    
    MolecularWeightsFactory mw_factory( ex, "Default", MolecularWeightsSettings{} );
    auto mw = mw_factory.get_instance();

    mw.modify_weights(lb);
    functional_type func = functional_type( ExchCXX::Backend::builtin, ExchCXX::Functional::PBE0, ExchCXX::Spin::Unpolarized );
    XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", 
        integrator_kernel, lwd_kernel, reduction_kernel );
    auto integrator = integrator_factory.get_instance( func, lb );

    OneDFTSettings onedft_settings;
    onedft_settings.model = onedft_model_path;

    auto [ EXC, VXC, VXCz ] = integrator.eval_exc_vxc_onedft( P, Pz, onedft_settings );
    auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
    auto VXCz_diff_nrm = ( VXCz - VXCz_ref ).norm();
    CHECK( EXC == Approx( EXC_ref ) );
    CHECK( VXC_diff_nrm / basis.nbf() < 1e-7 );
    CHECK( VXCz_diff_nrm / basis.nbf() < 1e-10 );
    // Check if the integrator propagates state correctly
    {
    auto [ EXC1, VXC1, VXCz1 ] = integrator.eval_exc_vxc_onedft( P, Pz, onedft_settings );
    CHECK( EXC1 == Approx( EXC_ref ) );
    auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
    auto VXCz1_diff_nrm = ( VXCz1 - VXCz_ref ).norm();
    CHECK( VXC1_diff_nrm / basis.nbf() < 1e-7 ); // TODO: Check this
    CHECK( VXCz1_diff_nrm / basis.nbf() < 1e-10 );
    }
}

void test_integrator(std::string reference_file, std::string onedft_model_path, bool use_cpu = true, bool use_gpu = true) {

#ifdef GAUXC_HAS_DEVICE
    auto rt = DeviceRuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9);
#else
    auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
#endif

#ifdef GAUXC_HAS_HOST
    if (use_cpu) {
        SECTION( "Host" ) {
        test_onedft_integrator( ExecutionSpace::Host, rt,
            reference_file, onedft_model_path );
        }
    }
#endif

#ifdef GAUXC_HAS_DEVICE
    if (use_gpu) {
        SECTION( "Device" ) {
        SECTION( "Incore - MPI Reduction" ) {
            test_onedft_integrator( ExecutionSpace::Device, rt,
                reference_file, onedft_model_path );
        }
        #ifdef GAUXC_HAS_CUTLASS
        SECTION( "Incore - MPI Reduction - CUTLASS" ) {
            test_onedft_integrator( ExecutionSpace::Device, rt,
            reference_file, onedft_model_path,
            "Default", "Default", "Scheme1-CUTLASS" );
        }
        #endif
        #ifdef GAUXC_HAS_NCCL
        SECTION( "Incore - NCCL Reduction" ) {
            test_onedft_integrator( ExecutionSpace::Device, rt,
            reference_file, onedft_model_path,
            "Default", "NCCL" );
        }
        #endif
        }
    }
#endif
}
    
TEST_CASE( "OneDFT", "[onedft]" ) {
    SECTION( " HE / def2-qzvp / tpss.fun" ) {
        test_integrator( GAUXC_REF_DATA_PATH "/onedft_he_def2qzvp_tpss_uks.hdf5", GAUXC_ONEDFT_MODEL_PATH "/tpss.fun" );
        }
    SECTION( " HE / def2-qzvp / pbe.fun" ) {
        test_integrator( GAUXC_REF_DATA_PATH "/onedft_he_def2qzvp_pbe_uks.hdf5", GAUXC_ONEDFT_MODEL_PATH "/pbe.fun" );
        }
    SECTION( " HE / def2-qzvp / lda.fun" ) {
        test_integrator( GAUXC_REF_DATA_PATH "/onedft_he_def2qzvp_lda_uks.hdf5", GAUXC_ONEDFT_MODEL_PATH "/lda.fun" );
        }
}

  // Regression guard for communicator-consistency bugs in OneDFT MPI paths.
  // Why this exists:
  // 1) The default MPI test runs with 2 ranks and RuntimeEnvironment(MPI_COMM_WORLD),
  //    which can mask bugs where collectives accidentally use MPI_COMM_WORLD.
  // 2) The real failure mode appears when RuntimeEnvironment is built from a
  //    subgroup communicator (comm_size < world_size), e.g. NGROUPS-style splits.
  // 3) This test forces that topology and exercises OneDFT so gather/scatter and
  //    metadata sizing must all agree on rt.comm().
TEST_CASE( "OneDFT MPI Subgroup", "[onedft][mpi][subcomm]" ) {
#ifdef GAUXC_HAS_MPI
  int world_rank = 0;
  int world_size = 1;
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );

  if( world_size < 3 ) {
    SUCCEED( "Requires at least 3 MPI ranks" );
    return;
  }

  MPI_Comm subcomm = MPI_COMM_NULL;
  const int color = world_rank % 2;
  MPI_Comm_split( MPI_COMM_WORLD, color, world_rank, &subcomm );

  int sub_size = 0;
  MPI_Comm_size( subcomm, &sub_size );
  CHECK( sub_size < world_size );

#ifdef GAUXC_HAS_DEVICE
  auto rt = DeviceRuntimeEnvironment( GAUXC_MPI_CODE(subcomm,) 0.9 );
#else
  auto rt = RuntimeEnvironment( GAUXC_MPI_CODE(subcomm) );
#endif

#ifdef GAUXC_HAS_HOST
  SECTION( "Host / HE / def2-qzvp / tpss.fun" ) {
    test_onedft_integrator( ExecutionSpace::Host, rt,
      GAUXC_REF_DATA_PATH "/onedft_he_def2qzvp_tpss_uks.hdf5",
      GAUXC_ONEDFT_MODEL_PATH "/tpss.fun" );
  }
#endif

#ifdef GAUXC_HAS_DEVICE
  SECTION( "Device / HE / def2-qzvp / tpss.fun" ) {
    test_onedft_integrator( ExecutionSpace::Device, rt,
      GAUXC_REF_DATA_PATH "/onedft_he_def2qzvp_tpss_uks.hdf5",
      GAUXC_ONEDFT_MODEL_PATH "/tpss.fun" );
  }
#endif

  MPI_Comm_free( &subcomm );
#else
  SUCCEED( "MPI disabled" );
#endif
}

#if defined(GAUXC_HAS_HOST) && defined(GAUXC_HAS_DEVICE)
void test_onedft_grad_host_device( std::string reference_file,
    std::string onedft_model_path ) {

    using matrix_type = Eigen::MatrixXd;
    Molecule mol;
    BasisSet<double> basis;
    read_hdf5_record( mol,   reference_file, "/MOLECULE" );
    read_hdf5_record( basis, reference_file, "/BASIS"    );

    HighFive::File file( reference_file, HighFive::File::ReadOnly );
    auto dset = file.getDataSet( "/DENSITY" );
    auto dims = dset.getDimensions();
    matrix_type P( dims[0], dims[1] );
    dset.read( P.data() );

    // We can only call OneDFT with UKS (two channels) for now, so we create a dummy Pz matrix to 
    // satisfy the interface.
    matrix_type Ps = P;
    matrix_type Pz = matrix_type::Zero( dims[0], dims[1] );

    auto mg = MolGridFactory::create_default_molgrid( mol, PruningScheme::Unpruned,
        BatchSize(512), RadialQuad::MuraKnowles, AtomicGridSizeDefault::UltraFineGrid );

    functional_type func = functional_type( ExchCXX::Backend::builtin,
        ExchCXX::Functional::PBE0, ExchCXX::Spin::Unpolarized );

    OneDFTSettings onedft_settings;
    onedft_settings.model = onedft_model_path;

#ifdef GAUXC_HAS_DEVICE
    auto rt = DeviceRuntimeEnvironment( GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9 );
#else
    auto rt = RuntimeEnvironment( GAUXC_MPI_CODE(MPI_COMM_WORLD) );
#endif

    auto eval_grad = [&]( ExecutionSpace ex ) {
        LoadBalancerFactory lb_factory( ex, "Default" );
        auto lb = lb_factory.get_instance( rt, mol, mg, basis );
        MolecularWeightsFactory mw_factory( ex, "Default", MolecularWeightsSettings{} );
        auto mw = mw_factory.get_instance();
        mw.modify_weights( lb );
        XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated",
            "Default", "Default", "Default" );
        auto integrator = integrator_factory.get_instance( func, lb );
        return integrator.eval_exc_grad_onedft( Ps, Pz, onedft_settings );
    };

    auto host_grad = eval_grad( ExecutionSpace::Host );
    auto dev_grad  = eval_grad( ExecutionSpace::Device );

    const size_t n = 3 * mol.size();
    REQUIRE( host_grad.size() == n );
    REQUIRE( dev_grad.size()  == n );

    // We want a non-zero gradient to make sure we have something to compare against.
    double host_squared_norm = 0.0;
    for( size_t i = 0; i < n; ++i ) host_squared_norm += host_grad[i] * host_grad[i];
    CHECK( std::sqrt( host_squared_norm ) > 1e-3 );

    // The device gradient must match the host gradient component-wise.
    for( size_t i = 0; i < n; ++i ) {
        CHECK( dev_grad[i] == Approx( host_grad[i] ).margin( 1e-6 ) );
    }
}

TEST_CASE( "OneDFT EXC Gradient", "[onedft][grad]" ) {
    SECTION( " H2O2 / def2-tzvp / tpss.fun" ) {
        test_onedft_grad_host_device(
            GAUXC_REF_DATA_PATH "/h2o2_def2-tzvp.hdf5",
            GAUXC_ONEDFT_MODEL_PATH "/tpss.fun" );
    }
}
#endif

#include <gauxc/xc_integrator/replicated/impl.hpp>
// Include the OneDFT utility header for reorder helper functions
#include "../../src/xc_integrator/integrator_util/onedft_util.hpp"

TEST_CASE( "Atom Reorder Permutation", "[onedft][reorder]" ) {

  // Scenario: 2 ranks, 3 atoms
  // Rank 0 has: atom0=2pts, atom1=3pts, atom2=1pt  (6 pts total)
  // Rank 1 has: atom0=1pt,  atom1=0pts, atom2=2pts (3 pts total)
  //
  // Rank-ordered layout (what MPI_Gatherv produces):
  //   [r0_a0(2), r0_a1(3), r0_a2(1), r1_a0(1), r1_a1(0), r1_a2(2)]
  //   indices: 0 1 | 2 3 4 | 5 | 6 | | 7 8
  //
  // Atom-ordered layout (what we want):
  //   [a0_r0(2), a0_r1(1), a1_r0(3), a1_r1(0), a2_r0(1), a2_r1(2)]
  //   indices: 0 1 | 2 | 3 4 5 | | 6 | 7 8

  int natoms = 3;
  int world_size = 2;
  // all_rank_atom_sizes: [rank0_atom0, rank0_atom1, rank0_atom2, rank1_atom0, rank1_atom1, rank1_atom2]
  std::vector<int64_t> all_rank_atom_sizes = {2, 3, 1, 1, 0, 2};
  std::vector<int> sendcounts = {6, 3};
  std::vector<int> displs = {0, 6};

  SECTION("Permutation correctness") {
    auto [perm, inv_perm] = GauXC::build_atom_reorder_perm(
      all_rank_atom_sizes, sendcounts, displs, natoms, world_size);

    REQUIRE(perm.size() == 9);
    REQUIRE(inv_perm.size() == 9);

    // Expected mapping:
    // rank-ordered idx -> atom-ordered idx
    // r0_a0: src 0->dst 0, src 1->dst 1
    // r0_a1: src 2->dst 3, src 3->dst 4, src 4->dst 5
    // r0_a2: src 5->dst 6
    // r1_a0: src 6->dst 2
    // r1_a1: (empty)
    // r1_a2: src 7->dst 7, src 8->dst 8
    CHECK(perm[0] == 0);
    CHECK(perm[1] == 1);
    CHECK(perm[2] == 3);
    CHECK(perm[3] == 4);
    CHECK(perm[4] == 5);
    CHECK(perm[5] == 6);
    CHECK(perm[6] == 2);
    CHECK(perm[7] == 7);
    CHECK(perm[8] == 8);

    // Round-trip: inv_perm[perm[i]] == i
    for (int64_t i = 0; i < 9; ++i) {
      CHECK(inv_perm[perm[i]] == i);
    }
  }

  SECTION("Strided permutation with stride 3 (coords)") {
    auto [perm, inv_perm] = GauXC::build_atom_reorder_perm(
      all_rank_atom_sizes, sendcounts, displs, natoms, world_size);

    // 9 points, stride 3 -> 27 doubles
    // Rank-ordered data: point i has values [i*10, i*10+1, i*10+2]
    std::vector<double> src(27);
    for (int i = 0; i < 9; ++i) {
      src[i*3]   = i * 10.0;
      src[i*3+1] = i * 10.0 + 1.0;
      src[i*3+2] = i * 10.0 + 2.0;
    }

    std::vector<double> dst(27, -1.0);
    GauXC::apply_strided_permutation(src.data(), dst.data(), perm, 9, 3);

    // Verify: dst[perm[i]*3..] should equal src[i*3..]
    for (int i = 0; i < 9; ++i) {
      int64_t j = perm[i];
      CHECK(dst[j*3]   == src[i*3]);
      CHECK(dst[j*3+1] == src[i*3+1]);
      CHECK(dst[j*3+2] == src[i*3+2]);
    }
  }

  SECTION("Round-trip: forward then inverse restores original") {
    auto [perm, inv_perm] = GauXC::build_atom_reorder_perm(
      all_rank_atom_sizes, sendcounts, displs, natoms, world_size);

    // Stride 1 (like grid_weights)
    std::vector<double> original(9);
    for (int i = 0; i < 9; ++i) original[i] = i * 1.5 + 0.7;

    // Forward: rank-ordered -> atom-ordered
    std::vector<double> atom_ordered(9);
    GauXC::apply_strided_permutation(original.data(), atom_ordered.data(), perm, 9, 1);

    // Inverse: atom-ordered -> rank-ordered
    std::vector<double> restored(9);
    GauXC::apply_strided_permutation(atom_ordered.data(), restored.data(), inv_perm, 9, 1);

    for (int i = 0; i < 9; ++i) {
      CHECK(restored[i] == Approx(original[i]));
    }
  }

  SECTION("Single rank is identity permutation") {
    // With 1 rank, no reorder is needed
    std::vector<int64_t> single_rank_sizes = {2, 3, 1};
    std::vector<int> sc = {6};
    std::vector<int> dp = {0};

    auto [perm, inv_perm] = GauXC::build_atom_reorder_perm(
      single_rank_sizes, sc, dp, 3, 1);

    for (int64_t i = 0; i < 6; ++i) {
      CHECK(perm[i] == i);
      CHECK(inv_perm[i] == i);
    }
  }

  SECTION("reorder_to_atom_order / reorder_to_rank_order round-trip") {
    auto [perm, inv_perm] = GauXC::build_atom_reorder_perm(
      all_rank_atom_sizes, sendcounts, displs, natoms, world_size);
    int64_t npts = 9;

    // Create synthetic interleaved data matching the real layouts
    std::vector<double> weights(npts), den(npts*2), coords(npts*3), dden(npts*6), tau_v(npts*2);
    for (int64_t i = 0; i < npts; ++i) {
      weights[i] = i * 0.1;
      den[i*2] = i * 1.0; den[i*2+1] = i * 1.0 + 100;
      coords[i*3] = i; coords[i*3+1] = i+0.1; coords[i*3+2] = i+0.2;
      for (int c = 0; c < 6; ++c) dden[i*6+c] = i * 10.0 + c;
      tau_v[i*2] = i * 5.0; tau_v[i*2+1] = i * 5.0 + 50;
    }
    // Save originals
    auto orig_weights = weights, orig_den = den, orig_coords = coords;
    auto orig_dden = dden, orig_tau = tau_v;

    // Forward: rank-order → atom-order
    GauXC::reorder_to_atom_order(weights, den, coords, dden, tau_v, perm, npts);

    // Verify data actually changed (perm is non-trivial)
    bool any_different = false;
    for (int64_t i = 0; i < npts && !any_different; ++i)
      if (weights[i] != orig_weights[i]) any_different = true;
    CHECK(any_different);

    // Now simulate the gradient path: convert interleaved atom-ordered data
    // to channel-first layout (as mpi_scatter_onedft_outputs does)
    std::vector<double> grad_den(npts*2), grad_dden(npts*6), grad_tau(npts*2);
    // Channel-first: [alpha(npts) | beta(npts)]
    for (int64_t i = 0; i < npts; ++i) {
      grad_den[i] = den[i*2];              // alpha channel
      grad_den[npts+i] = den[i*2+1];       // beta channel
      grad_tau[i] = tau_v[i*2];
      grad_tau[npts+i] = tau_v[i*2+1];
      // dden channel-first: [dXa(npts)|dYa|dZa|dXb|dYb|dZb]
      for (int c = 0; c < 6; ++c)
        grad_dden[c*npts+i] = dden[i*6+c];
    }

    // Inverse: atom-order → rank-order
    GauXC::reorder_to_rank_order(grad_den, grad_dden, grad_tau, inv_perm, npts, true, true);

    // Verify round-trip: channel-first rank-ordered should match original interleaved
    for (int64_t i = 0; i < npts; ++i) {
      CHECK(grad_den[i] == Approx(orig_den[i*2]));
      CHECK(grad_den[npts+i] == Approx(orig_den[i*2+1]));
      CHECK(grad_tau[i] == Approx(orig_tau[i*2]));
      CHECK(grad_tau[npts+i] == Approx(orig_tau[i*2+1]));
      for (int c = 0; c < 6; ++c)
        CHECK(grad_dden[c*npts+i] == Approx(orig_dden[i*6+c]));
    }
  }
}