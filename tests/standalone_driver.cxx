#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/runtime_environment.hpp>
#include <gauxc/molecular_weights.hpp>
#include <gauxc/molgrid/defaults.hpp>

#include <gauxc/external/hdf5.hpp>
#include <highfive/H5File.hpp>
#include "ini_input.hpp"
#include <gauxc/exceptions.hpp>
#define EIGEN_DONT_VECTORIZE
#define EIGEN_NO_CUDA
#include <Eigen/Core>

using namespace GauXC;
using namespace ExchCXX;

int main(int argc, char** argv) {

#ifdef GAUXC_ENABLE_MPI
  MPI_Init( NULL, NULL );
#endif
  {

    // Set up runtimes
    #ifdef GAUXC_ENABLE_DEVICE
    DeviceRuntimeEnvironment rt( GAUXC_MPI_CODE(MPI_COMM_WORLD), 0.9 );
    #else
    RuntimeEnvironment rt(GAUXC_MPI_CODE(MPI_COMM_WORLD));
    #endif
    auto world_rank = rt.comm_rank();
    auto world_size = rt.comm_size();

    std::vector< std::string > opts( argc );
    for( int i = 0; i < argc; ++i ) opts[i] = argv[i];

    //std::string ref_file  = opts.at(1);
    auto input_file = opts.at(1);
    INIFile input(input_file);

    // Require Ref file
    auto ref_file = input.getData<std::string>("GAUXC.REF_FILE");

    // Optional Args
    std::string grid_spec = "ULTRAFINE";
    std::string prune_spec = "UNPRUNED";
    std::string lb_exec_space_str = "Host";
    std::string int_exec_space_str = "Host";
    std::string integrator_kernel = "Default";
    std::string lwd_kernel         = "Default";
    std::string reduction_kernel   = "Default";
    double      basis_tol = 1e-10;
    std::string func_spec = "PBE0";
    bool integrate_den      = false;
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

    if( input.containsData("GAUXC.PRUNING_SCHEME") ) {
      prune_spec = input.getData<std::string>("GAUXC.PRUNING_SCHEME");
      string_to_upper( prune_spec );
    }

    if( input.containsData("GAUXC.BASIS_TOL") )
      basis_tol = input.getData<double>("GAUXC.BASIS_TOL");

    if( input.containsData("GAUXC.FUNC") ) {
      func_spec = input.getData<std::string>("GAUXC.FUNC");
      string_to_upper( func_spec );
    }

    if( input.containsData("GAUXC.INTEGRATE_DEN" ) )
      integrate_den = input.getData<bool>("GAUXC.INTEGRATE_DEN");
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

    #ifdef GAUXC_ENABLE_DEVICE
    std::map< std::string, ExecutionSpace > exec_space_map = {
      { "HOST",   ExecutionSpace::Host },
      { "DEVICE", ExecutionSpace::Device }
    };

    auto lb_exec_space = exec_space_map.at(lb_exec_space_str);
    auto int_exec_space = exec_space_map.at(int_exec_space_str);
    #else
    auto lb_exec_space  = ExecutionSpace::Host;
    auto int_exec_space = ExecutionSpace::Host;
    #endif

    if( !world_rank ) {
      std::cout << "DRIVER SETTINGS: " << std::endl
                << "  REF_FILE          = " << ref_file << std::endl
                << "  GRID              = " << grid_spec << std::endl
                << "  PRUNING SCHEME    = " << prune_spec << std::endl
                << "  BASIS_TOL         = " << basis_tol << std::endl
                << "  FUNCTIONAL        = " << func_spec << std::endl
                << "  LB_EXEC_SPACE     = " << lb_exec_space_str << std::endl
                << "  INT_EXEC_SPACE    = " << int_exec_space_str << std::endl
                << "  INTEGRATOR_KERNEL = " << integrator_kernel << std::endl
                << "  LWD_KERNEL        = " << lwd_kernel << std::endl
                << "  REDUCTION_KERNEL  = " << reduction_kernel << std::endl
                << "  VXC (?)           = " << std::boolalpha << integrate_vxc << std::endl
                << "  EXX (?)           = " << std::boolalpha << integrate_exx << std::endl
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
    //  std::cout << x.Z.get() << ", " 
    //    << x.x << ", " << x.y << ", " << x.z << std::endl;
    //}

    // Construct MolGrid / MolMeta
    std::map< std::string, AtomicGridSizeDefault > mg_map = {
      {"FINE",      AtomicGridSizeDefault::FineGrid},
      {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
      {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid}
    };

    std::map< std::string, PruningScheme > prune_map = {
      {"UNPRUNED", PruningScheme::Unpruned},
      {"ROBUST",   PruningScheme::Robust},
      {"TREUTLER", PruningScheme::Treutler}
    };

    auto mg = MolGridFactory::create_default_molgrid(mol, 
     prune_map.at(prune_spec), BatchSize(4096), 
     RadialQuad::MuraKnowles, mg_map.at(grid_spec));

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
    auto lb = lb_factory.get_shared_instance( rt, mol, mg, basis);

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

    // Apply molecular partition weights
    MolecularWeightsFactory mw_factory( int_exec_space, "Default", MolecularWeightsSettings{} );
    auto mw = mw_factory.get_instance();
    mw.modify_weights(*lb);

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
    std::vector<double> EXC_GRAD_ref(3*mol.size());
    size_t N_EL_ref = MolMeta(mol).sum_atomic_charges();
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
        try {
          dset = file.getDataSet("/VXC");
          dset.read( VXC_ref.data() );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "Could Not Find Reference VXC" << std::endl;
          }
          VXC_ref.fill(0);
        }

        try {
          dset = file.getDataSet("/EXC");
          dset.read( &EXC_ref );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "Could Not Find Reference EXC" << std::endl;
          }
          EXC_ref = 0.;
        }
      }

      if( integrate_exc_grad ) {
        try {
          dset = file.getDataSet("EXC_GRAD");
          auto xc_grad_dims = dset.getDimensions();
          if( xc_grad_dims[0] != mol.size() or xc_grad_dims[1] != 3 )
            throw std::runtime_error("Incorrect dims for EXC_GRAD");
          dset.read( EXC_GRAD_ref.data() );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "Could Not Find Reference EXC_GRAD" << std::endl;
          }
          std::fill( EXC_GRAD_ref.begin(), EXC_GRAD_ref.end(), 0. );
        }
      }

      if( integrate_exx ) {
        try {
          dset = file.getDataSet("/K");
          dset.read( K_ref.data() );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "Could Not Find Reference K" << std::endl;
          }
          K_ref.fill(0);
        }
      }

    }

    //std::cout << "NBF = " << basis.nbf() << std::endl;
    
#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    auto xc_int_start = std::chrono::high_resolution_clock::now();

    matrix_type VXC, K;
    double EXC, N_EL;

    if( integrate_den ) {
      N_EL = integrator.integrate_den( P );
      std::cout << std::scientific << std::setprecision(12);
      std::cout << "N_EL = " << N_EL << std::endl;
    } else {
      N_EL = N_EL_ref;
    }

    if( integrate_vxc ) {
      std::tie(EXC, VXC) = integrator.eval_exc_vxc( P );
      std::cout << std::scientific << std::setprecision(12);
      std::cout << "EXC = " << EXC << std::endl;
    } else {
      EXC = EXC_ref;
      VXC = VXC_ref;
    }

    std::vector<double> EXC_GRAD;
    if( integrate_exc_grad ) {
      EXC_GRAD = integrator.eval_exc_grad( P );
      std::cout << "EXC Gradient:" << std::endl;
      std::cout << std::scientific << std::setprecision(6);
      for( auto iAt = 0; iAt < mol.size(); ++iAt ) {
        std::cout << "  " 
                  << std::setw(16) << EXC_GRAD[3*iAt + 0] 
                  << std::setw(16) << EXC_GRAD[3*iAt + 1] 
                  << std::setw(16) << EXC_GRAD[3*iAt + 2] 
                  << std::endl;
      }
    }

    if( integrate_exx ) {
      K = integrator.eval_exx(P, sn_link_settings);
      matrix_type K_tmp = 0.5 * (K + K.transpose());
      K = K_tmp;
    } else                K = K_ref;

    //std::cout << (K).block(0,0,5,5) << std::endl << std::endl;
    //std::cout << (K_ref).block(0,0,5,5) << std::endl << std::endl;
    //std::cout << (K - K_ref).block(0,0,5,5) << std::endl << std::endl;
#ifdef GAUXC_ENABLE_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    auto xc_int_end   = std::chrono::high_resolution_clock::now();
    double xc_int_dur = std::chrono::duration<double>( xc_int_end - xc_int_start ).count();

#ifdef GAUXC_ENABLE_MPI
    util::MPITimer mpi_lb_timings( MPI_COMM_WORLD, lb->get_timings() );
    util::MPITimer mpi_xc_timings( MPI_COMM_WORLD, integrator.get_timings() );
    util::MPITimer mpi_weight_timings( MPI_COMM_WORLD, mw.get_timings() );
#endif
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

      std::cout << "MolecularWeights Timings" << std::endl;
      for( const auto& [name, dur] : mw.get_timings().all_timings() ) {
        #ifdef GAUXC_ENABLE_MPI
        const auto avg     = mpi_weight_timings.get_avg_duration(name).count();
        const auto min     = mpi_weight_timings.get_min_duration(name).count();
        const auto max     = mpi_weight_timings.get_max_duration(name).count();
        const auto std_dev = mpi_weight_timings.get_std_dev(name).count();
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

      if( integrate_den ) {
      std::cout << "N_EL (ref)        = " << (double)N_EL_ref << std::endl;
      std::cout << "N_EL (calc)       = " << N_EL     << std::endl;
      std::cout << "N_EL Diff         = " << std::abs(N_EL_ref - N_EL) / N_EL_ref 
                                         << std::endl;
      }

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

      if(integrate_exc_grad) {
      double exc_grad_ref_nrm(0.), exc_grad_calc_nrm(0.), exc_grad_diff_nrm(0.);
      for( auto i = 0; i < 3*mol.size(); ++i ) {
        const auto ref_val = EXC_GRAD_ref[i];
        const auto clc_val = EXC_GRAD[i];
        const auto dif_val = std::abs(ref_val - clc_val);
        exc_grad_ref_nrm  += ref_val*ref_val;
        exc_grad_calc_nrm += clc_val*clc_val;
        exc_grad_diff_nrm += dif_val*dif_val;
      }

      exc_grad_ref_nrm  = std::sqrt(exc_grad_ref_nrm);
      exc_grad_calc_nrm = std::sqrt(exc_grad_calc_nrm);
      exc_grad_diff_nrm = std::sqrt(exc_grad_diff_nrm);
      std::cout << "| EXC_GRAD (ref)  | = " << exc_grad_ref_nrm << std::endl;
      std::cout << "| EXC_GRAD (calc) | = " << exc_grad_calc_nrm << std::endl;
      std::cout << "| EXC_GRAD (diff) | = " << exc_grad_diff_nrm << std::endl;
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
        HighFive::DataSpace grad_space( mol.size(), 3 );
        dset = file.createDataSet<double>( "/EXC_GRAD", grad_space );
        dset.write_raw( EXC_GRAD.data() );
      }
    }

  }
#ifdef GAUXC_ENABLE_MPI
  MPI_Finalize();
#endif

}
