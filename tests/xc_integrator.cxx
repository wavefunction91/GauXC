#include "ut_common.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/oop_xc_integrator/impl.hpp>
#include <gauxc/oop_xc_integrator/integrator_factory.hpp>

using namespace GauXC;

#ifdef GAUXC_ENABLE_MPI
  #define GAUXC_MPI_ARG   MPI_Comm comm,
  #define GAUXC_MPI_PARAM comm,
#else
  #define GAUXC_MPI_ARG                 
  #define GAUXC_MPI_PARAM      
#endif


void test_xc_integrator( ExecutionSpace ex, GAUXC_MPI_ARG Molecule mol, std::string integrator_kernel = "Default",  const bool check_state_propagation = true ) 
{

  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Default");
  auto lb = lb_factory.get_instance(GAUXC_MPI_PARAM mol, mg, basis);

  functional_type func( ExchCXX::Backend::builtin, ExchCXX::Functional::PBE0, ExchCXX::Spin::Unpolarized );

  using matrix_type = Eigen::MatrixXd;
  //auto integrator = make_default_integrator<matrix_type>( ex, comm, func, basis, lb );
  XCIntegratorFactory<matrix_type> integrator_factory( ex, "Replicated", integrator_kernel, "Default" );
  auto integrator = integrator_factory.get_instance( func, lb );


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


  if( ex == ExecutionSpace::Host ) {
    std::cout << "EXC = " << EXC << std::endl;
    auto EXC_GRAD = integrator.eval_exc_grad( P );
    for( auto x : EXC_GRAD ) std::cout << x << std::endl;
  }
}


TEST_CASE( "Benzene / PBE0 / cc-pVDZ", "[xc-integrator]" ) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  Molecule mol  = make_benzene();

#ifdef GAUXC_ENABLE_HOST
  SECTION( "Host" ) {
    test_xc_integrator( ExecutionSpace::Host, GAUXC_MPI_PARAM mol );
  }
#endif

#ifdef GAUXC_ENABLE_CUDA
  SECTION( "Device" ) {
    SECTION( "Default" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_PARAM mol, "Default" );
    }
    SECTION( "ShellBatched" ) {
      test_xc_integrator( ExecutionSpace::Device, GAUXC_MPI_PARAM mol, "ShellBatched" );
    }
  }
#endif

}
