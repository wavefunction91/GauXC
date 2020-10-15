#include "standards.hpp"
#include "basis/parse_basis.hpp"

#include <random>
#include <algorithm>

#include <fstream>
#include <gauxc/external/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include "eigen3_matrix_serialization.hpp"

#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>

#include <Eigen/Core>

using namespace GauXC;
using namespace ExchCXX;

int main(int argc, char** argv) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Init( NULL, NULL );
  int world_rank, world_size;
  MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &world_size );
#else
  int world_rank = 0;
  int world_size = 1;
#endif
  {

    std::vector< std::string > opts( argc );
    for( int i = 0; i < argc; ++i ) opts[i] = argv[i];

    std::string test_case = opts.at(1);
    std::string basis_set = opts.at(2);
    std::string ref_file  = opts.at(3);

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

    //std::cout << "Basis" << std::endl;
    //for( auto i = 0; i < basis.size(); ++i ) {
    //  const auto& sh = basis[i];
    //  std::cout << "CEN = " << sh.O()[0] << ", " << sh.O()[1] << ", " << sh.O()[2] << std::endl;
    //  std::cout << "L = " << sh.l() << std::endl;
    //  std::cout << "CR = " << sh.cutoff_radius() << std::endl;
    //  std::cout << "PRIMS" << std::endl;
    //  for( auto p = 0; p < sh.nprim(); ++p )
    //    std::cout << "  " << sh.alpha()[p] << ", " << sh.coeff()[p] << std::endl;
    //  std::cout << std::endl;
    //}

    // Setup load balancer
#ifdef GAUXC_ENABLE_MPI
    auto lb = std::make_shared<LoadBalancer>(MPI_COMM_WORLD, mol, mg, basis, meta);
#else
    auto lb = std::make_shared<LoadBalancer>(mol, mg, basis, meta);
#endif

    // Setup XC functional
    functional_type func( Backend::builtin, XCFunctional::Functional::PBE0, Spin::Unpolarized );

    // Setup Integrator
    using matrix_type = Eigen::MatrixXd;
#ifdef GAUXC_ENABLE_MPI
    XCIntegrator<matrix_type> integrator( ExecutionSpace::Device, MPI_COMM_WORLD, func, 
                                          basis, lb );
#else
    XCIntegrator<matrix_type> integrator( ExecutionSpace::Device, func, basis, lb );
#endif

    matrix_type P,VXC_ref;
    double EXC_ref;
    {
      std::ifstream infile( ref_file, std::ios::binary );

      if( !infile.good() ) throw std::runtime_error(ref_file + " not found");
      cereal::BinaryInputArchive ar(infile);
      ar( EXC_ref, P, VXC_ref );
    }

    //std::cout << "NBF = " << basis.nbf() << std::endl;
    
#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    auto xc_int_start = std::chrono::high_resolution_clock::now();

    auto [ EXC, VXC ] = integrator.eval_exc_vxc( P );

#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    auto xc_int_end   = std::chrono::high_resolution_clock::now();
    double xc_int_dur = std::chrono::duration<double>( xc_int_end - xc_int_start ).count();

    if( !world_rank ) {

      std::cout << "Load Balancer Timings" << std::endl;
      for( const auto& [name, dur] : lb->get_timings().all_timings() ) {
        std::cout << "  " << std::setw(30) << name << ": " 
                  << std::setw(10) << dur.count() << " ms" << std::endl;
      }

      std::cout << "Integrator Timings" << std::endl;
      for( const auto& [name, dur] : integrator.get_timings().all_timings() ) {
        std::cout << "  " << std::setw(30) << name << ": " 
                  << std::setw(10) << dur.count() << " ms" << std::endl;
      }

      std::cout << std::scientific << std::setprecision(16);

      std::cout << "XC Int Duration  = " << xc_int_dur << " s" << std::endl;

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
#ifdef GAUXC_ENABLE_MPI
  MPI_Finalize();
#endif

}
