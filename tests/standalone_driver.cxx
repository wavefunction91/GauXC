#include "standards.hpp"
#include <random>
#include <algorithm>

#include <fstream>
#include <cereal/archives/binary.hpp>
#include "eigen3_matrix_serialization.hpp"

#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>

#include <Eigen/Core>

using namespace GauXC;
using namespace ExchCXX;

int main(int argc, char** argv) {

  MPI_Init( NULL, NULL );
  int world_rank, world_size;
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
  {

    //std::string test_case = "benzene";
    //std::string basis_set = "cc-pvdz";
    std::string test_case = "taxol";
    std::string basis_set = "6-31gd";
    std::transform( test_case.begin(), test_case.end(), test_case.begin(),
                    [](const auto c){ return std::tolower(c); } );
    std::transform( basis_set.begin(), basis_set.end(), basis_set.begin(),
                    [](const auto c){ return std::tolower(c); } );


    // Construct Molecule
    Molecule mol;
    if( test_case.find("benzene") != std::string::npos )
      mol = make_benzene();
    else if( test_case.find("taxol") != std::string::npos )
      mol = make_taxol();
    else if( test_case.find("ubiquitin") != std::string::npos )
      mol = make_ubiquitin();
    else
      throw std::runtime_error("Unknown Test Case");

    MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);
    auto meta = std::make_shared<MolMeta>( mol );


    // Construct BasisSet
    BasisSet<double> basis; 
    if( basis_set.find("6-31gd") != std::string::npos ) 
      basis = std::move(make_631Gd( mol, SphericalType(false) ));
    else if( basis_set.find("cc-pvdz") != std::string::npos ) 
      basis = std::move(make_ccpvdz( mol, SphericalType(true) ));
    else
      throw std::runtime_error("Unknown Basis Set");

    for( auto& sh : basis ) 
      sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
    basis.generate_shell_to_ao();

    std::cout << "Basis" << std::endl;
    for( auto i = 0; i < basis.size(); ++i ) {
      const auto& sh = basis[i];
      std::cout << "CEN = " << sh.O()[0] << ", " << sh.O()[1] << ", " << sh.O()[2] << std::endl;
      std::cout << "L = " << sh.l() << std::endl;
      std::cout << "CR = " << sh.cutoff_radius() << std::endl;
      std::cout << "PRIMS" << std::endl;
      for( auto p = 0; p < sh.nprim(); ++p )
        std::cout << "  " << sh.alpha()[p] << ", " << sh.coeff()[p] << std::endl;
      std::cout << std::endl;
    }

    // Setup load balancer
    auto lb = std::make_shared<LoadBalancer>(MPI_COMM_WORLD, mol, mg, basis, meta);

    // Setup XC functional
    functional_type func( XCFunctional::Functional::PBE0, Spin::Unpolarized );

    // Setup Integrator
    using matrix_type = Eigen::MatrixXd;
    XCIntegrator<matrix_type> integrator( ExecutionSpace::Device, MPI_COMM_WORLD, func, 
                                          basis, lb );

    matrix_type P,VXC_ref;
    double EXC_ref;
    {
      std::string ref_file = test_case + "_pbe0_" + basis_set + "_ufg_ssf.bin";
      std::ifstream infile( ref_file, std::ios::binary );

      if( !infile.good() ) throw std::runtime_error(ref_file + " not found");
      cereal::BinaryInputArchive ar(infile);
      ar( EXC_ref, P, VXC_ref );
    }

    std::cout << "NBF = " << basis.nbf() << std::endl;
    auto [ EXC, VXC ] = integrator.eval_exc_vxc( P );

    if( !world_rank ) {
      std::cout << std::scientific << std::setprecision(16);
      std::cout << "EXC (ref)        = " << EXC_ref << std::endl;
      std::cout << "EXC (calc)       = " << EXC     << std::endl;
      std::cout << "EXC Diff         = " << std::abs(EXC_ref - EXC) / EXC_ref 
                                         << std::endl;

      std::cout << "| VXC (ref)  |_F = " << VXC_ref.norm() << std::endl;
      std::cout << "| VXC (calc) |_F = " << VXC.norm() << std::endl;
      std::cout << "RMS VXC Diff     = " << (VXC_ref - VXC).norm() / basis.nbf()
                                         << std::endl;
    }
  
  }
  MPI_Finalize();

}
