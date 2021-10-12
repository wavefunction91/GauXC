#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/util/div_ceil.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include "ini_input.hpp"
#include <gauxc/exceptions.hpp>
#define EIGEN_DONT_VECTORIZE
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

    //std::string ref_file  = opts.at(1);
    auto input_file = opts.at(1);
    INIFile input(input_file);

    // Require Ref file
    auto ref_file = input.getData<std::string>("GAUXC.REF_FILE");

    // Optional Args
    std::string grid_spec = "ULTRAFINE";
    std::string lb_exec_space_str = "Host";
    std::string int_exec_space_str = "Host";
    std::string integrator_kernel = "Default";
    std::string lwd_kernel         = "Default";
    std::string reduction_kernel   = "Default";
    double      basis_tol = 1e-10;
    std::string func_spec = "PBE0";
    bool integrate_vxc      = true;
    bool integrate_exx      = false;
    bool integrate_exc_grad = false;

    auto string_to_upper = []( auto& str ) {
      std::transform( str.begin(), str.end(), str.begin(), ::toupper );
    };

    if( input.containsData("GAUXC.GRID") ) {
      grid_spec = input.getData<std::string>("GAUXC.GRID");
      string_to_upper( grid_spec );
    }

    if( input.containsData("GAUXC.BASIS_TOL") )
      basis_tol = input.getData<double>("GAUXC.BASIS_TOL");

    if( input.containsData("GAUXC.FUNC") ) {
      func_spec = input.getData<std::string>("GAUXC.FUNC");
      string_to_upper( func_spec );
    }

    if( input.containsData("GAUXC.INTEGRATE_VXC" ) )
      integrate_vxc = input.getData<bool>("GAUXC.INTEGRATE_VXC");
    if( input.containsData("GAUXC.INTEGRATE_EXX" ) )
      integrate_exx = input.getData<bool>("GAUXC.INTEGRATE_EXX");
    if( input.containsData("GAUXC.INTEGRATE_EXC_GRAD" ) )
      integrate_exc_grad = input.getData<bool>("GAUXC.INTEGRATE_EXC_GRAD");

    if( input.containsData("GAUXC.LB_EXEC_SPACE") )
      lb_exec_space_str = input.getData<std::string>("GAUXC.LB_EXEC_SPACE");
    if( input.containsData("GAUXC.INT_EXEC_SPACE") )
      int_exec_space_str = input.getData<std::string>("GAUXC.INT_EXEC_SPACE");

    if( input.containsData("GAUXC.INTEGRATOR_KERNEL") )
      integrator_kernel = input.getData<std::string>("GAUXC.INTEGRATOR_KERNEL");
    if( input.containsData("GAUXC.LWD_KERNEL") )
      lwd_kernel = input.getData<std::string>("GAUXC.LWD_KERNEL");
    if( input.containsData("GAUXC.REDUCTION_KERNEL") )
      reduction_kernel = input.getData<std::string>("GAUXC.REDUCTION_KERNEL");

    string_to_upper( lb_exec_space_str );
    string_to_upper( int_exec_space_str );
    string_to_upper( integrator_kernel );
    string_to_upper( lwd_kernel );
    string_to_upper( reduction_kernel );

    std::map< std::string, ExecutionSpace > exec_space_map = {
      { "HOST",   ExecutionSpace::Host },
      #ifdef GAUXC_ENABLE_DEVICE
      { "DEVICE", ExecutionSpace::Device }
      #endif
    };

    auto lb_exec_space = exec_space_map.at(lb_exec_space_str);
    auto int_exec_space = exec_space_map.at(int_exec_space_str);

    if( !world_rank ) {
      std::cout << "DRIVER SETTINGS: " << std::endl
                << "  REF_FILE       = " << ref_file << std::endl
                << "  GRID           = " << grid_spec << std::endl
                << "  BASIS_TOL      = " << basis_tol << std::endl
                << "  FUNCTIONAL     = " << func_spec << std::endl
                << "  LB_EXEC_SPACE  = " << lb_exec_space_str << std::endl
                << "  INT_EXEC_SPACE = " << int_exec_space_str << std::endl
                << "  VXC (?)        = " << std::boolalpha << integrate_vxc << std::endl
                << "  EXX (?)        = " << std::boolalpha << integrate_exx << std::endl
                << std::endl;
    }

    IntegratorSettingsSNLinK sn_link_settings;
    if( input.containsData("EXX.TOL_E") )
      sn_link_settings.energy_tol = input.getData<double>("EXX.TOL_E");
    if( input.containsData("EXX.TOL_K") )
      sn_link_settings.k_tol = input.getData<double>("EXX.TOL_K");

    // Read Molecule
    Molecule mol;
    read_hdf5_record( mol, ref_file, "/MOLECULE" );
    //std::cout << "Molecule" << std::endl;
    //for( auto x : mol ) {
    //  std::cout << x.Z.get() << ", " << x.x << ", " << x.y << ", " << x.z << std::endl;
    //}

    // Construct MolGrid / MolMeta
    std::map< std::string, AtomicGridSizeDefault > mg_map = {
      {"FINE",      AtomicGridSizeDefault::FineGrid},
      {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
      {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid}
    };

    MolGrid mg( mg_map.at(grid_spec), mol);

    // Read BasisSet
    BasisSet<double> basis; 
    read_hdf5_record( basis, ref_file, "/BASIS" );
    //std::cout << basis.nbf() << std::endl;

    for( auto& sh : basis ){ 
      sh.set_shell_tolerance( basis_tol );
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
    //LoadBalancerFactory lb_factory(ExecutionSpace::Device, "Default");
    //LoadBalancerFactory lb_factory(ExecutionSpace::Host, "Replicated-FillIn");
    LoadBalancerFactory lb_factory( lb_exec_space, "Replicated");
    auto lb = lb_factory.get_shared_instance( GAUXC_MPI_CODE(MPI_COMM_WORLD,) mol, 
      mg, basis);

    if(0){
      auto& tasks = lb->get_tasks();
      size_t total_npts = 0, total_npts_fixed = 0;
      for( auto&& task : tasks ) {
        auto npts = task.points.size();
        auto npts_fix = util::div_ceil( npts, 32 ) * 32;
        total_npts += npts;
        total_npts_fixed += npts_fix;
      }
      std::cout << "TOTAL NPTS = " << total_npts << " , FIXED = " << total_npts_fixed << std::endl;
    }

    // Setup XC functional
    functional_type func( Backend::builtin, functional_map.value(func_spec), 
      Spin::Unpolarized );

    // Setup Integrator
    using matrix_type = Eigen::MatrixXd;
    XCIntegratorFactory<matrix_type> integrator_factory( int_exec_space , 
      "Replicated", integrator_kernel, lwd_kernel, reduction_kernel );
    auto integrator = integrator_factory.get_instance( func, lb );

    // Read in reference data
    matrix_type P,VXC_ref,K_ref;
    double EXC_ref;
    {
      HighFive::File file( ref_file, HighFive::File::ReadOnly );
      auto dset = file.getDataSet("/DENSITY");
      auto dims = dset.getDimensions();
      P       = matrix_type( dims[0], dims[1] );
      VXC_ref = matrix_type( dims[0], dims[1] );
      K_ref   = matrix_type( dims[0], dims[1] );

      if( P.rows() != P.cols() ) throw std::runtime_error("Density Must Be Square");
      if( P.rows() != basis.nbf() ) 
        throw std::runtime_error("Density Not Compatible With Basis");

      dset.read( P.data() );

      if( integrate_vxc ) {
        dset = file.getDataSet("/VXC");
        dset.read( VXC_ref.data() );

        dset = file.getDataSet("/EXC");
        dset.read( &EXC_ref );
      }

      if( integrate_exx ) {
        dset = file.getDataSet("/K");
        dset.read( K_ref.data() );
      }

    }

    //std::cout << "NBF = " << basis.nbf() << std::endl;
    
#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    auto xc_int_start = std::chrono::high_resolution_clock::now();

    matrix_type VXC, K;
    double EXC;

      std::cout << std::scientific << std::setprecision(6);
    if( integrate_vxc ) {
      std::tie(EXC, VXC) = integrator.eval_exc_vxc( P );
      std::cout << "EXC = " << EXC << std::endl;
    } else {
      EXC = EXC_ref;
      VXC = VXC_ref;
    }

    std::vector<double> EXC_GRAD;
    if( integrate_exc_grad ) {
      EXC_GRAD = integrator.eval_exc_grad( P );
      std::cout << "EXC Gradient:" << std::endl;
      for( auto iAt = 0; iAt < mol.size(); ++iAt ) {
        std::cout << "  " 
                  << std::setw(16) << EXC_GRAD[3*iAt + 0] 
                  << std::setw(16) << EXC_GRAD[3*iAt + 1] 
                  << std::setw(16) << EXC_GRAD[3*iAt + 2] 
                  << std::endl;
      }
    }

    if( integrate_exx ) K = integrator.eval_exx(P, sn_link_settings);
    else                K = K_ref;

    //std::cout << (K).block(0,0,5,5) << std::endl << std::endl;
    //std::cout << (K_ref).block(0,0,5,5) << std::endl << std::endl;
    //std::cout << (K - K_ref).block(0,0,5,5) << std::endl << std::endl;
#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    auto xc_int_end   = std::chrono::high_resolution_clock::now();
    double xc_int_dur = std::chrono::duration<double>( xc_int_end - xc_int_start ).count();

    util::MPITimer mpi_lb_timings( MPI_COMM_WORLD, lb->get_timings() );
    util::MPITimer mpi_xc_timings( MPI_COMM_WORLD, integrator.get_timings() );
    if( !world_rank ) {

      std::cout << std::scientific << std::setprecision(5) << std::endl;
      std::cout << "Load Balancer Timings" << std::endl;
      for( const auto& [name, dur] : lb->get_timings().all_timings() ) {
        #ifdef GAUXC_ENABLE_MPI
        const auto avg     = mpi_lb_timings.get_avg_duration(name).count();
        const auto min     = mpi_lb_timings.get_min_duration(name).count();
        const auto max     = mpi_lb_timings.get_max_duration(name).count();
        const auto std_dev = mpi_lb_timings.get_std_dev(name).count();
        #endif
        std::cout << "  " << std::setw(30) << name << ": " 
        #ifdef GAUXC_ENABLE_MPI
                  << "AVG = " << std::setw(12) << avg << " ms, " 
                  << "MIN = " << std::setw(12) << min << " ms, " 
                  << "MAX = " << std::setw(12) << max << " ms, " 
                  << "STDDEV = " << std::setw(12) << std_dev << " ms" << std::endl;
        #else
                  << std::setw(12) << dur.count() << " ms" << std::endl;
        #endif
      }

      std::cout << "Integrator Timings" << std::endl;
      for( const auto& [name, dur] : integrator.get_timings().all_timings() ) {
        #ifdef GAUXC_ENABLE_MPI
        const auto avg     = mpi_xc_timings.get_avg_duration(name).count();
        const auto min     = mpi_xc_timings.get_min_duration(name).count();
        const auto max     = mpi_xc_timings.get_max_duration(name).count();
        const auto std_dev = mpi_xc_timings.get_std_dev(name).count();
        #endif
        std::cout << "  " << std::setw(30) << name << ": " 
        #ifdef GAUXC_ENABLE_MPI
                  << "AVG = " << std::setw(12) << avg << " ms, " 
                  << "MIN = " << std::setw(12) << min << " ms, " 
                  << "MAX = " << std::setw(12) << max << " ms, " 
                  << "STDDEV = " << std::setw(12) << std_dev << " ms" << std::endl;
        #else
                  << std::setw(12) << dur.count() << " ms" << std::endl;
        #endif
      }

      std::cout << std::scientific << std::setprecision(14);

      std::cout << "XC Int Duration  = " << xc_int_dur << " s" << std::endl;

      if( integrate_vxc ) {
      std::cout << "EXC (ref)        = " << EXC_ref << std::endl;
      std::cout << "EXC (calc)       = " << EXC     << std::endl;
      std::cout << "EXC Diff         = " << std::abs(EXC_ref - EXC) / EXC_ref 
                                         << std::endl;

      std::cout << "| VXC (ref)  |_F = " << VXC_ref.norm() << std::endl;
      std::cout << "| VXC (calc) |_F = " << VXC.norm() << std::endl;
      std::cout << "RMS VXC Diff     = " << (VXC_ref - VXC).norm() / basis.nbf()
                                         << std::endl;
      }

      if( integrate_exx ) {
      std::cout << "| K (ref)  |_F = " << K_ref.norm() << std::endl;
      std::cout << "| K (calc) |_F = " << K.norm() << std::endl;
      std::cout << "RMS K Diff     = " << (K_ref - K).norm() / basis.nbf()
                                         << std::endl;
      }
    }

    // Dump out new file
    if( input.containsData("GAUXC.OUTFILE") ) {
      // Create File
      auto outfname = input.getData<std::string>("GAUXC.OUTFILE");
      { HighFive::File( outfname, HighFive::File::Truncate ); }

      // Write molecule
      write_hdf5_record( mol, outfname, "/MOLECULE" );

      // Write Basis
      write_hdf5_record( basis, outfname, "/BASIS" );

      // Write out matrices
      HighFive::File file( outfname, HighFive::File::ReadWrite );
      HighFive::DataSpace mat_space( basis.nbf(), basis.nbf() );
      HighFive::DataSpace sca_space( 1 );

      auto dset = file.createDataSet<double>( "/DENSITY", mat_space );
      dset.write_raw( P.data() );

      if( integrate_vxc ) {
        dset = file.createDataSet<double>( "/VXC", mat_space );
        dset.write_raw( VXC.data() );

        dset = file.createDataSet<double>( "/EXC", sca_space );
        dset.write_raw( &EXC );
      }

      if( integrate_exx ) {
        dset = file.createDataSet<double>( "/K", mat_space );
        dset.write_raw( K.data() );
      }

      if( integrate_exc_grad ) {
        HighFive::DataSpace grad_space( 3*mol.size() );
        dset = file.createDataSet<double>( "/EXC_GRAD", grad_space );
        dset.write_raw( EXC_GRAD.data() );
      }
    }

  }
#ifdef GAUXC_ENABLE_MPI
  MPI_Finalize();
#endif

}
