#include "ut_common.hpp"
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>

#include <Eigen/Core>

using namespace GauXC;

TEST_CASE( "Benzene / PBE0 / cc-pVDZ", "[xc-integrator]" ) {

  MPI_Comm comm          = MPI_COMM_WORLD;
  Molecule mol           = make_benzene();
  BasisSet<double> basis = make_ccpvdz( mol, SphericalType(true) );

  for( auto& sh : basis ) 
    sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
  basis.generate_shell_to_ao();

  MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);

  auto meta = std::make_shared<MolMeta>( mol );

  auto lb = std::make_shared<LoadBalancer>(comm, mol, mg, basis, meta);

  functional_type func( ExchCXX::XCFunctional::Functional::PBE0, ExchCXX::Spin::Unpolarized );

  using matrix_type = Eigen::MatrixXd;
  XCIntegrator<matrix_type> integrator( comm, func, basis, lb );


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

  auto VXC_diff_nrm = ( VXC - VXC_ref ).norm();
  CHECK( VXC_diff_nrm / basis.nbf() < 1e-10 ); 
}
