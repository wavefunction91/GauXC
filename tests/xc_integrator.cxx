#include "ut_common.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/new_xc_integrator/impl.hpp>
#include <gauxc/new_xc_integrator/integrator_factory.hpp>

using namespace GauXC;



#ifdef GAUXC_ENABLE_MPI
void test_xc_integrator( ExecutionSpace ex, MPI_Comm comm, Molecule mol, const bool check_state_propagation = true ) 
#else
void test_xc_integrator( ExecutionSpace ex, Molecule mol, const bool check_state_propagation = true ) 
#endif
{

  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  auto meta = std::make_shared<MolMeta>( mol );

#ifdef GAUXC_ENABLE_MPI
  auto lb = std::make_shared<LoadBalancer>(comm, mol, mg, basis, meta);
#else
  auto lb = std::make_shared<LoadBalancer>(mol, mg, basis, meta);
#endif

  functional_type func( ExchCXX::Backend::builtin, ExchCXX::Functional::PBE0, ExchCXX::Spin::Unpolarized );

  using matrix_type = Eigen::MatrixXd;
#ifdef GAUXC_ENABLE_MPI
  auto integrator = make_default_integrator<matrix_type>( ex, comm, func, basis, lb );
#else
  auto integrator = make_default_integartor<matrix_type>( ex, func, basis, lb );
#endif


  matrix_type P,VXC_ref;
  double EXC_ref;

  {
    std::string ref_file = GAUXC_REF_DATA_PATH "/benzene_pbe0_cc-pvdz_ufg_ssf.bin";
    std::ifstream infile( ref_file, std::ios::binary );
    cereal::BinaryInputArchive ar(infile);
    ar( EXC_ref, P, VXC_ref );
  }

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
  if( check_state_propagation ) {
    auto [ EXC1, VXC1 ] = integrator.eval_exc_vxc( P );
    CHECK( EXC1 == Approx( EXC_ref ) );
    auto VXC1_diff_nrm = ( VXC1 - VXC_ref ).norm();
    CHECK( VXC1_diff_nrm / basis.nbf() < 1e-10 ); 
  }
}


TEST_CASE( "Benzene / PBE0 / cc-pVDZ", "[xc-integrator]" ) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  Molecule mol  = make_benzene();

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host" ) {
#ifdef GAUXC_ENABLE_MPI
    test_xc_integrator( ExecutionSpace::Host, comm, mol );
#else
    test_xc_integrator( ExecutionSpace::Host, mol );
#endif
  }
#endif

#ifdef GAUXC_ENABLE_CUDA
  SECTION( "Device" ) {
#ifdef GAUXC_ENABLE_MPI
    test_xc_integrator( ExecutionSpace::Device, comm, mol );
#else
    test_xc_integrator( ExecutionSpace::Device, mol );
#endif
  }
#endif

}
