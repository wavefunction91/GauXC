#include <gauxc/xc_integrator.hpp>
#include <gauxc/oop_xc_integrator/impl.hpp>
#include <gauxc/oop_xc_integrator/integrator_factory.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
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

    std::string ref_file  = opts.at(1);

    // Read Molecule
    Molecule mol;
    read_hdf5_record( mol, ref_file, "/MOLECULE" );

    // Construct MolGrid / MolMeta
    MolGrid mg(AtomicGridSizeDefault::UltraFineGrid, mol);
    auto meta = std::make_shared<MolMeta>( mol );

    // Read BasisSet
    BasisSet<double> basis; 
    read_hdf5_record( basis, ref_file, "/BASIS" );

    for( auto& sh : basis ){ 
      //sh.set_shell_tolerance( std::numeric_limits<double>::epsilon() );
      sh.set_shell_tolerance( 1e-10 );
    }

    //std::cout << "Basis" << std::endl;
    //for( auto& sh : basis ) std::cout << sh << std::endl;

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
    LoadBalancerFactory lb_factory(ExecutionSpace::Device, "Default");
    auto lb = lb_factory.get_shared_instance( GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, mg, basis);

    // Setup XC functional
    functional_type func( Backend::builtin, Functional::PBE0, Spin::Unpolarized );

    // Setup Integrator
    using matrix_type = Eigen::MatrixXd;
    XCIntegratorFactory<matrix_type> integrator_factory( ExecutionSpace::Device, 
      "Replicated", "Default", "Default" );
    auto integrator = integrator_factory.get_instance( func, lb );

    // Read in reference data
    matrix_type P,VXC_ref;
    double EXC_ref;
    {
      HighFive::File file( ref_file, HighFive::File::ReadOnly );
      auto dset = file.getDataSet("/DENSITY");
      auto dims = dset.getDimensions();
      P       = matrix_type( dims[0], dims[1] );
      VXC_ref = matrix_type( dims[0], dims[1] );

      dset.read( P.data() );
      dset = file.getDataSet("/VXC");
      dset.read( VXC_ref.data() );

      dset = file.getDataSet("/EXC");
      dset.read( &EXC_ref );
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
