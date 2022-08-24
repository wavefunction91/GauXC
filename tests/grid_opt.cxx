#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/basisset_map.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include "ini_input.hpp"
#include <gauxc/exceptions.hpp>
#define EIGEN_DONT_VECTORIZE
#include <Eigen/Core>

#include <integratorxx/quadratures/muraknowles.hpp>
#include <integratorxx/quadratures/treutleraldrichs.hpp>
#include <integratorxx/quadratures/lebedev_laikov.hpp>
#include <integratorxx/composite_quadratures/pruned_spherical_quadrature.hpp>
#include <gauxc/grid_factory.hpp>
#include <gauxc/molgrid/defaults.hpp>

#include <chrono>

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

    // Get input file
    auto input_file = opts.at(1);
    INIFile input(input_file);

    // Require Ref file
    auto ref_file = input.getData<std::string>("GAUXC.REF_FILE");
    //auto ref_ne   = input.getData<double>("GAUXC.REF_NE");

    // Read Molecule
    Molecule mol;
    read_hdf5_record( mol, ref_file, "/MOLECULE" );

    double ref_ne = 0.0;
    for( auto& atom : mol ) {
      ref_ne += atom.Z.get();
    }

    // Read BasisSet
    BasisSet<double> basis; 
    read_hdf5_record( basis, ref_file, "/BASIS" );

    for( auto& sh : basis ){ sh.set_shell_tolerance( 1e-10 ); }


    // Setup load balancer factory
    LoadBalancerFactory lb_factory( ExecutionSpace::Host, "Replicated");

    // Setup Integrator factory
    using matrix_type = Eigen::MatrixXd;
    XCIntegratorFactory<matrix_type> integrator_factory( ExecutionSpace::Host, 
      "Replicated", "Default", "Default", "Default" );

    // Setup Dummy XC functional
    functional_type func;

    // Read in reference density
    matrix_type P;
    {
      HighFive::File file( ref_file, HighFive::File::ReadOnly );
      auto dset = file.getDataSet("/DENSITY");
      auto dims = dset.getDimensions();
      P       = matrix_type( dims[0], dims[1] );

      if( P.rows() != P.cols() ) 
        throw std::runtime_error("Density Must Be Square");
      if( P.rows() != basis.nbf() ) 
        throw std::runtime_error("Density Not Compatible With Basis");

      dset.read( P.data() );
    }


#if 1
    UnprunedAtomicGridSpecification cno_unp_spec {
      RadialQuad::MuraKnowles,
      RadialSize(100),
      RadialScale(7.0),
      AngularSize(974)
    };

    UnprunedAtomicGridSpecification h_unp_spec {
      RadialQuad::MuraKnowles,
      RadialSize(100),
      RadialScale(5.0),
      AngularSize(974)
    };
    auto c_unp_grid = AtomicGridFactory::generate_grid(cno_unp_spec);
    auto n_unp_grid = AtomicGridFactory::generate_grid(cno_unp_spec);
    auto o_unp_grid = AtomicGridFactory::generate_grid(cno_unp_spec);
    auto h_unp_grid = AtomicGridFactory::generate_grid(h_unp_spec);

    atomic_grid_map unp_molmap = {
      { AtomicNumber(1), h_unp_grid },
      { AtomicNumber(6), c_unp_grid },
      { AtomicNumber(7), n_unp_grid },
      { AtomicNumber(8), o_unp_grid }
    };
#else
    auto unp_molmap = MolGridFactory
#endif

    std::cout << "Unpruned" << std::endl;
    // Unpruned Integration
    {
      auto st = std::chrono::high_resolution_clock::now();
      MolGrid mg(unp_molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      auto en = std::chrono::high_resolution_clock::now();
      std::cout << std::scientific << std::setprecision(16);
      const auto err = std::abs(N_EL-ref_ne);
      std::cout << "NE = " << N_EL << ", " << err << ", " << err/ref_ne << std::endl;
      std::cout << std::chrono::duration<double>(en-st).count() << std::endl;
    
    }

   // Setup pruning
   #if 0
   #if 0
   std::vector<PruningRegion> pruning_regions = {
     {0ul, 25ul,   AngularSize(170)},
     {25ul, 50ul,  AngularSize(266)},
     {50ul, 100ul, AngularSize(974)}
   };

   PrunedAtomicGridSpecification cno_pru_spec {
     RadialQuad::MuraKnowles,
     RadialSize(100),
     RadialScale(7.0),
     pruning_regions
   };

   PrunedAtomicGridSpecification h_pru_spec {
     RadialQuad::MuraKnowles,
     RadialSize(100),
     RadialScale(5.0),
     pruning_regions
   };
   #else
   auto cno_pru_spec = create_pruned_spec( PruningScheme::Robust, cno_unp_spec );
   auto h_pru_spec   = create_pruned_spec( PruningScheme::Robust, h_unp_spec   );
   #endif

    auto c_pru_grid = AtomicGridFactory::generate_grid(cno_pru_spec);
    auto n_pru_grid = AtomicGridFactory::generate_grid(cno_pru_spec);
    auto o_pru_grid = AtomicGridFactory::generate_grid(cno_pru_spec);
    auto h_pru_grid = AtomicGridFactory::generate_grid(h_pru_spec);

    atomic_grid_map pru_molmap = {
      { AtomicNumber(1), h_pru_grid },
      { AtomicNumber(6), c_pru_grid },
      { AtomicNumber(7), n_pru_grid },
      { AtomicNumber(8), o_pru_grid }
    };
    #else
    #if 0
    atomic_grid_map pru_molmap;
    for( auto& atom : mol ) {
      if(!pru_molmap.count(atom.Z)) {
        pru_molmap.emplace(atom.Z,AtomicGridFactory::generate_grid(
          MolGridFactory::create_default_pruned_grid_spec(
            PruningScheme::Robust, atom.Z, RadialQuad::MuraKnowles,
            RadialSize(100), AngularSize(974)
          )
        ));
      }
    }
    #else
    auto pru_molmap = MolGridFactory::create_default_gridmap(
      mol, PruningScheme::Robust, RadialQuad::MuraKnowles,
      RadialSize(100), AngularSize(974) );
    #endif
    #endif

    std::cout << "Pruned" << std::endl;
    // Pruned Integration
    {
      auto st = std::chrono::high_resolution_clock::now();
      MolGrid mg(pru_molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      auto en = std::chrono::high_resolution_clock::now();
      std::cout << std::scientific << std::setprecision(16);
      const auto err = std::abs(N_EL-ref_ne);
      std::cout << "NE = " << N_EL << ", " << err << ", " << err/ref_ne << std::endl;
      std::cout << std::chrono::duration<double>(en-st).count() << std::endl;
    
    }


  }
#ifdef GAUXC_ENABLE_MPI
  MPI_Finalize();
#endif

}
