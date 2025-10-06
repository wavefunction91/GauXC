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