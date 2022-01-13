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
    auto ref_ne   = input.getData<double>("GAUXC.REF_NE");

    // Read Molecule
    Molecule mol;
    read_hdf5_record( mol, ref_file, "/MOLECULE" );

    // Read BasisSet
    BasisSet<double> basis; 
    read_hdf5_record( basis, ref_file, "/BASIS" );

    for( auto& sh : basis ){ sh.set_shell_tolerance( 1e-14 ); }


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
    IntegratorXX::MuraKnowles<double,double> c_rad( 100, 7.0 );
    IntegratorXX::MuraKnowles<double,double> n_rad( 100, 7.0 );
    IntegratorXX::MuraKnowles<double,double> o_rad( 100, 7.0 );
    IntegratorXX::MuraKnowles<double,double> h_rad( 100, 5.0 );
    IntegratorXX::LebedevLaikov<double> ang(974);
    IntegratorXX::LebedevLaikov<double> ang_m6(266);
    //IntegratorXX::LebedevLaikov<double> ang(3470);
    //IntegratorXX::LebedevLaikov<double> ang_m6(1454);
    IntegratorXX::LebedevLaikov<double> ang_7(170);

#if 0
    using sph_type = IntegratorXX::SphericalQuadrature<
      IntegratorXX::MuraKnowles<double,double>,
      IntegratorXX::LebedevLaikov<double>
    >;

    Grid c_unp_grid(std::make_shared<sph_type>( c_rad, ang ));
    Grid n_unp_grid(std::make_shared<sph_type>( n_rad, ang ));
    Grid o_unp_grid(std::make_shared<sph_type>( o_rad, ang ));
    Grid h_unp_grid(std::make_shared<sph_type>( h_rad, ang ));
#else
    auto c_unp_grid = AtomicGridFactory::generate_unpruned_grid( c_rad, ang );
    auto n_unp_grid = AtomicGridFactory::generate_unpruned_grid( n_rad, ang );
    auto o_unp_grid = AtomicGridFactory::generate_unpruned_grid( o_rad, ang );
    auto h_unp_grid = AtomicGridFactory::generate_unpruned_grid( h_rad, ang );
#endif

#if 0
    IntegratorXX::RadialGridPartition 
      rgp_c( c_rad, 0, ang_7, 26, ang_m6, 51, ang );
    IntegratorXX::RadialGridPartition 
      rgp_n( n_rad, 0, ang_7, 26, ang_m6, 51, ang );
    IntegratorXX::RadialGridPartition 
      rgp_o( o_rad, 0, ang_7, 26, ang_m6, 51, ang );
    IntegratorXX::RadialGridPartition 
      rgp_h( h_rad, 0, ang_7, 26, ang_m6, 51, ang );

    using pruned_sph_type = IntegratorXX::PrunedSphericalQuadrature<
      IntegratorXX::MuraKnowles<double,double>,
      IntegratorXX::LebedevLaikov<double>
    >;

    Grid c_pru_grid(std::make_shared<pruned_sph_type>( c_rad, rgp_c ));
    Grid n_pru_grid(std::make_shared<pruned_sph_type>( n_rad, rgp_c ));
    Grid o_pru_grid(std::make_shared<pruned_sph_type>( o_rad, rgp_c ));
    Grid h_pru_grid(std::make_shared<pruned_sph_type>( h_rad, rgp_h ));
#else
    auto c_pru_grid = AtomicGridFactory::generate_pruned_grid( c_rad, ang, ang_m6, ang_7 );
    auto n_pru_grid = AtomicGridFactory::generate_pruned_grid( n_rad, ang, ang_m6, ang_7 );
    auto o_pru_grid = AtomicGridFactory::generate_pruned_grid( o_rad, ang, ang_m6, ang_7 );
    auto h_pru_grid = AtomicGridFactory::generate_pruned_grid( h_rad, ang, ang_m6, ang_7 );
#endif

    std::cout << "Unpruned" << std::endl;
    // Unpruned Integration
    {
      auto st = std::chrono::high_resolution_clock::now();
      atomic_grid_map molmap = {
        { AtomicNumber(1), h_unp_grid },
        { AtomicNumber(6), c_unp_grid },
        { AtomicNumber(7), n_unp_grid },
        { AtomicNumber(8), o_unp_grid }
      };
      MolGrid mg(molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      auto en = std::chrono::high_resolution_clock::now();
      std::cout << std::scientific << std::setprecision(4);
      const auto err = std::abs(N_EL-ref_ne);
      std::cout << "NE = " << err << ", " << err/ref_ne << std::endl;
      std::cout << std::chrono::duration<double>(en-st).count() << std::endl;
    
    }

    std::cout << "Pruned" << std::endl;
    // Pruned Integration
    {
      auto st = std::chrono::high_resolution_clock::now();
      atomic_grid_map molmap = {
        { AtomicNumber(1), h_pru_grid },
        { AtomicNumber(6), c_pru_grid },
        { AtomicNumber(7), n_pru_grid },
        { AtomicNumber(8), o_pru_grid }
      };
      MolGrid mg(molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      auto en = std::chrono::high_resolution_clock::now();
      std::cout << std::scientific << std::setprecision(4);
      const auto err = std::abs(N_EL-ref_ne);
      std::cout << "NE = " << err << ", " << err/ref_ne << std::endl;
      std::cout << std::chrono::duration<double>(en-st).count() << std::endl;
    
    }
#else

    BasisSetMap bs_map( basis, mol );

    const size_t nr = 125;
    const size_t na = 3470;

#if 1
    const auto rad_quad = RadialQuad::MuraKnowles;
    IntegratorXX::MuraKnowles<double,double> mk_grid(nr,1.0);
    auto max_raw_r = mk_grid.points().back();
    const double def_c_r = 7.0;
    const double def_h_r = 5.0;
#else
    const auto rad_quad = RadialQuad::TreutlerAldrichs;
    IntegratorXX::TreutlerAldrichs<double,double> mk_grid(nr,1.0);
    auto max_raw_r = mk_grid.points().back();
    const double def_c_r = 0.8;
    const double def_h_r = 1.1;
#endif


    double avg_c_r = 0.; size_t n_c = 0;
    double avg_h_r = 0.; size_t n_h = 0;
    for( auto ish = 0; ish < basis.size(); ++ish ) {
      if( bs_map.shell_to_center(ish) == 0 ) {
        avg_c_r += basis[ish].cutoff_radius();
        n_c++;
      }

      if( bs_map.shell_to_center(ish) == 6 ) {
        avg_h_r += basis[ish].cutoff_radius();
        n_h++;
      }
    }

    avg_c_r /= n_c;
    avg_h_r /= n_h;

    std::cout << std::endl;
    double Rc = avg_c_r / max_raw_r;
    double Rh = avg_h_r / max_raw_r;
    
    std::cout << "Default Params" << std::endl;
    {
      Grid c_grid( rad_quad, RadialSize(nr), AngularSize(na), 
                   RadialScale(def_c_r) );

      Grid h_grid( rad_quad, RadialSize(nr), AngularSize(na), 
                   RadialScale(def_h_r) );


      atomic_grid_map molmap = {
        { AtomicNumber(1), h_grid },
        { AtomicNumber(6), c_grid }
      };
      MolGrid mg(molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      std::cout << std::scientific << std::setprecision(4);
      const auto err = std::abs(N_EL-ref_ne);
      std::cout << "NE = " << err << ", " << err/ref_ne << std::endl;
    }

#if 1
    std::cout << "Opt Params" << std::endl;
    std::cout << "C R = " << Rc << std::endl;
    std::cout << "H R = " << Rh << std::endl;
    {
      Grid c_grid( rad_quad, RadialSize(nr), AngularSize(na), 
                   RadialScale(Rc) );

      Grid h_grid( rad_quad, RadialSize(nr), AngularSize(na), 
                   RadialScale(Rh) );



      atomic_grid_map molmap = {
        { AtomicNumber(1), h_grid },
        { AtomicNumber(6), c_grid }
      };
      MolGrid mg(molmap);

      auto lb = lb_factory.get_shared_instance( 
        GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);
      auto integrator = integrator_factory.get_instance( func, lb );

      double N_EL = integrator.integrate_den( P );
      std::cout << std::scientific << std::setprecision(4);
      std::cout << "NE = " << std::abs(N_EL - ref_ne)/ref_ne << std::endl;
    }
#endif
#endif


  }
#ifdef GAUXC_ENABLE_MPI
  MPI_Finalize();
#endif

}
