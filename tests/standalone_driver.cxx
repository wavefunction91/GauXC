/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
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

#ifdef GAUXC_HAS_MPI
  MPI_Init( NULL, NULL );
#endif
  {

    // Set up runtimes
    #ifdef GAUXC_HAS_DEVICE
    auto rt = DeviceRuntimeEnvironment( GAUXC_MPI_CODE(MPI_COMM_WORLD,) 0.9 );
    #else
    auto rt = RuntimeEnvironment(GAUXC_MPI_CODE(MPI_COMM_WORLD));
    #endif
    auto world_rank = rt.comm_rank();
    auto world_size = rt.comm_size();

    std::vector< std::string > opts( argc );
    for( int i = 0; i < argc; ++i ) opts[i] = argv[i];

    auto input_file = opts.at(1);
    INIFile input(input_file);

    // Require Ref file
    auto ref_file = input.getData<std::string>("GAUXC.REF_FILE");

    // Optional Args
    std::string grid_spec          = "ULTRAFINE";
    std::string prune_spec         = "UNPRUNED";
    std::string lb_exec_space_str  = "Host";
    std::string int_exec_space_str = "Host";
    std::string integrator_kernel  = "Default";
    std::string lwd_kernel         = "Default";
    std::string reduction_kernel   = "Default";

    size_t      batch_size = 512;
    double      basis_tol  = 1e-10;
    std::string func_spec  = "PBE0";

    bool integrate_den      = false;
    bool integrate_vxc      = true;
    bool integrate_exx      = false;
    bool integrate_exc_grad = false;

    auto string_to_upper = []( auto& str ) {
      std::transform( str.begin(), str.end(), str.begin(), ::toupper );
    };

    #define OPTIONAL_KEYWORD(NAME,VAR,TYPE) \
    if( input.containsData(NAME) ) {        \
        VAR = input.getData<TYPE>(NAME);    \
    }

    OPTIONAL_KEYWORD( "GAUXC.GRID",              grid_spec,          std::string );
    OPTIONAL_KEYWORD( "GAUXC.FUNC",              func_spec,          std::string );
    OPTIONAL_KEYWORD( "GAUXC.PRUNING_SCHEME",    prune_spec,         std::string );
    OPTIONAL_KEYWORD( "GAUXC.LB_EXEC_SPACE",     lb_exec_space_str,  std::string );
    OPTIONAL_KEYWORD( "GAUXC.INT_EXEC_SPACE",    int_exec_space_str, std::string );
    OPTIONAL_KEYWORD( "GAUXC.INTEGRATOR_KERNEL", integrator_kernel,  std::string );
    OPTIONAL_KEYWORD( "GAUXC.LWD_KERNEL",        lwd_kernel,         std::string );
    OPTIONAL_KEYWORD( "GAUXC.REDUCTION_KERNEL",  reduction_kernel,   std::string );
    string_to_upper( grid_spec          );
    string_to_upper( func_spec          );
    string_to_upper( prune_spec         );
    string_to_upper( lb_exec_space_str  );
    string_to_upper( int_exec_space_str );
    string_to_upper( integrator_kernel  );
    string_to_upper( lwd_kernel         );
    string_to_upper( reduction_kernel   );

    OPTIONAL_KEYWORD( "GAUXC.BATCH_SIZE",     batch_size, size_t );
    OPTIONAL_KEYWORD( "GAUXC.BASIS_TOL",      basis_tol,  double );

    OPTIONAL_KEYWORD( "GAUXC.INTEGRATE_DEN",      integrate_den,      bool );
    OPTIONAL_KEYWORD( "GAUXC.INTEGRATE_VXC",      integrate_vxc,      bool );
    OPTIONAL_KEYWORD( "GAUXC.INTEGRATE_EXX",      integrate_exx,      bool );
    OPTIONAL_KEYWORD( "GAUXC.INTEGRATE_EXC_GRAD", integrate_exc_grad, bool );

    IntegratorSettingsSNLinK sn_link_settings;
    OPTIONAL_KEYWORD( "EXX.TOL_E", sn_link_settings.energy_tol, double );
    OPTIONAL_KEYWORD( "EXX.TOL_K", sn_link_settings.k_tol,      double );


    #ifdef GAUXC_HAS_DEVICE
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
      std::cout << std::boolalpha;
      std::cout << "DRIVER SETTINGS: " << std::endl
                << "  REF_FILE          = " << ref_file << std::endl
                << "  GRID              = " << grid_spec << std::endl
                << "  PRUNING_SCHEME    = " << prune_spec << std::endl
                << "  BATCH_SIZE        = " << batch_size << std::endl
                << "  BASIS_TOL         = " << basis_tol << std::endl
                << "  FUNCTIONAL        = " << func_spec << std::endl
                << "  LB_EXEC_SPACE     = " << lb_exec_space_str << std::endl
                << "  INT_EXEC_SPACE    = " << int_exec_space_str << std::endl
                << "  INTEGRATOR_KERNEL = " << integrator_kernel << std::endl
                << "  LWD_KERNEL        = " << lwd_kernel << std::endl
                << "  REDUCTION_KERNEL  = " << reduction_kernel << std::endl
                << "  DEN (?)           = " << integrate_den << std::endl
                << "  VXC (?)           = " << integrate_vxc << std::endl
                << "  EXX (?)           = " << integrate_exx << std::endl
                << "  EXC_GRAD (?)      = " << integrate_exc_grad << std::endl;
                if(integrate_exx) {
                  std::cout << "  EXX.TOL_E         = " 
                            << sn_link_settings.energy_tol << std::endl
                            << "  EXX.TOL_K         = " 
                            << sn_link_settings.k_tol << std::endl;
                }
                std::cout << std::endl;
    }






    // Read Molecule
    Molecule mol;
    read_hdf5_record( mol, ref_file, "/MOLECULE" );

    // Construct MolGrid / MolMeta
    std::map< std::string, AtomicGridSizeDefault > mg_map = {
      {"FINE",      AtomicGridSizeDefault::FineGrid},
      {"ULTRAFINE", AtomicGridSizeDefault::UltraFineGrid},
      {"SUPERFINE", AtomicGridSizeDefault::SuperFineGrid},
      {"GM3",       AtomicGridSizeDefault::GM3},
      {"GM5",       AtomicGridSizeDefault::GM5}
    };

    std::map< std::string, PruningScheme > prune_map = {
      {"UNPRUNED", PruningScheme::Unpruned},
      {"ROBUST",   PruningScheme::Robust},
      {"TREUTLER", PruningScheme::Treutler}
    };

    auto mg = MolGridFactory::create_default_molgrid(mol, 
     prune_map.at(prune_spec), BatchSize(batch_size), 
     RadialQuad::MuraKnowles, mg_map.at(grid_spec));

    // Read BasisSet
    BasisSet<double> basis; 
    read_hdf5_record( basis, ref_file, "/BASIS" );

    for( auto& sh : basis ){ 
      sh.set_shell_tolerance( basis_tol );
    }

    // Setup load balancer
    LoadBalancerFactory lb_factory( lb_exec_space, "Replicated");
    auto lb = lb_factory.get_shared_instance( rt, mol, mg, basis);

    // Apply molecular partition weights
    MolecularWeightsFactory mw_factory( int_exec_space, "Default", 
      MolecularWeightsSettings{} );
    auto mw = mw_factory.get_instance();
    mw.modify_weights(*lb);


    using matrix_type = Eigen::MatrixXd;
    // Read in reference data
    matrix_type P, Pz, Py, Px, VXC_ref, VXCz_ref, VXCy_ref, VXCx_ref, K_ref;
    double EXC_ref;
    std::vector<double> EXC_GRAD_ref(3*mol.size());
    bool rks = true, uks = false, gks = false;
    size_t N_EL_ref = MolMeta(mol).sum_atomic_charges();
    {
      HighFive::File file( ref_file, HighFive::File::ReadOnly );
      std::string den_str = "/DENSITY";
      std::string vxc_str = "/VXC";

      if (file.exist("/DENSITY_Z") ) {
        rks = false;
        den_str = "/DENSITY_SCALAR";
        vxc_str = "/VXC_SCALAR";
        if (file.exist("/DENSITY_Y") and file.exist("/DENSITY_X")) {
          gks = true;
        } else {
          uks = true;
        }
      }


      auto dset = file.getDataSet(den_str);
      auto dims = dset.getDimensions();
      P       = matrix_type( dims[0], dims[1] );
      VXC_ref = matrix_type( dims[0], dims[1] );
      K_ref   = matrix_type( dims[0], dims[1] );
      if (not rks) {
        Pz        = matrix_type( dims[0], dims[1] );
        VXCz_ref  = matrix_type( dims[0], dims[1] );
      }
      if (gks) {
        Py        = matrix_type( dims[0], dims[1] );
        VXCy_ref  = matrix_type( dims[0], dims[1] );
        Px        = matrix_type( dims[0], dims[1] );
        VXCx_ref  = matrix_type( dims[0], dims[1] );
      }

      if( P.rows() != P.cols() ) GAUXC_GENERIC_EXCEPTION("Density Must Be Square");
      if( P.rows() != basis.nbf() ) 
        GAUXC_GENERIC_EXCEPTION("Density Not Compatible With Basis");

      dset.read( P.data() );

      if( not rks) {
        dset = file.getDataSet("/DENSITY_Z");
        dset.read( Pz.data() );
      }
      if(gks) {
        dset = file.getDataSet("/DENSITY_Y");
        dset.read( Py.data() );
        dset = file.getDataSet("/DENSITY_X");
        dset.read( Px.data() );
      }
        

      if( integrate_vxc ) {
        try {
          dset = file.getDataSet(vxc_str);
          dset.read( VXC_ref.data() );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "** Warning: Could Not Find Reference VXC" << std::endl;
          }
          VXC_ref.fill(0);
        }
        if( not rks ) {
          try {
            dset = file.getDataSet("/VXC_Z");
            dset.read( VXCz_ref.data() );
          } catch(...) {
            if(world_rank == 0) {
              std::cout << "** Warning: Could Not Find Reference VXCz" << std::endl;
            }
            VXCz_ref.fill(0);
          }
        }
        if( gks ) {
          try {
            dset = file.getDataSet("/VXC_Y");
            dset.read( VXCy_ref.data() );
            dset = file.getDataSet("/VXC_X");
            dset.read( VXCx_ref.data() );
          } catch(...) {
            if(world_rank == 0) {
              std::cout << "** Warning: Could Not Find Reference VXCy/x" << std::endl;
            }
            VXCy_ref.fill(0);
            VXCx_ref.fill(0);
          }
        }

        try {
          dset = file.getDataSet("/EXC");
          dset.read( &EXC_ref );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "** Warning: Could Not Find Reference EXC" << std::endl;
          }
          EXC_ref = 0.;
        }

      }

      if( integrate_exc_grad ) {
        try {
          dset = file.getDataSet("EXC_GRAD");
          auto xc_grad_dims = dset.getDimensions();
          if( xc_grad_dims[0] != mol.size() or xc_grad_dims[1] != 3 )
            GAUXC_GENERIC_EXCEPTION("Incorrect dims for EXC_GRAD");
          dset.read( EXC_GRAD_ref.data() );
        } catch(...) {
          if(world_rank == 0) {
            std::cout << "** Warning: Could Not Find Reference EXC_GRAD" 
                      << std::endl;
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
            std::cout << "** Warning: Could Not Find Reference K" << std::endl;
          }
          K_ref.fill(0);
        }
      }

    }
    // Setup XC functional
    auto polar = (uks or gks) ? Spin::Polarized : Spin::Unpolarized;
    functional_type func( Backend::builtin, functional_map.value(func_spec), 
      polar );
    // Setup Integrator
    XCIntegratorFactory<matrix_type> integrator_factory( int_exec_space , 
      "Replicated", integrator_kernel, lwd_kernel, reduction_kernel );
    auto integrator = integrator_factory.get_instance( func, lb );
    
#ifdef GAUXC_HAS_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif
    auto xc_int_start = std::chrono::high_resolution_clock::now();

    matrix_type VXC, VXCz, VXCy, VXCx, K;
    double EXC, N_EL;

    std::cout << std::scientific << std::setprecision(12);
    if( integrate_den ) {
      if( (uks or gks) and !world_rank ) {
        std::cout << "Warning: integrate_den will only integrate the scalar density!" << std::endl;
      }
      N_EL = integrator.integrate_den( P );
      if(!world_rank) std::cout << "N_EL = " << N_EL << std::endl;
    } else {
      N_EL = N_EL_ref;
    }

    if( integrate_vxc ) {
      if( rks ) {
        std::tie(EXC, VXC) = integrator.eval_exc_vxc( P );
      }
      else if ( uks ) {
        std::tie(EXC, VXC, VXCz) = integrator.eval_exc_vxc( P, Pz );
      }
      else if ( gks ) {
        std::tie(EXC, VXC, VXCz, VXCy, VXCx) = integrator.eval_exc_vxc( P, Pz, Py, Px );
      }
      std::cout << std::scientific << std::setprecision(12);
      if(!world_rank) std::cout << "EXC = " << EXC << std::endl;
    } else {
      EXC = EXC_ref;
      VXC = VXC_ref;
      if( not rks ) {
        VXCz = VXCz_ref;
        if( gks ) {
          VXCy = VXCy_ref;
          VXCx = VXCx_ref;
        }
      }
    }

    std::vector<double> EXC_GRAD;
    if( integrate_exc_grad ) {
      if( rks ) {
        EXC_GRAD = integrator.eval_exc_grad( P );
      }
      else if( uks ) {
        std::cout << "Warning: eval_exc_grad + UKS NYI!" << std::endl;
        //EXC_GRAD = integrator.eval_exc_grad( P, Pz );
      }
      else if( gks ) {
        std::cout << "Warning: eval_exc_grad + GKS NYI!" << std::endl;
        //EXC_GRAD = integrator.eval_exc_grad( P, Pz, Py, Px );
      }
      if(!world_rank) {
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
    }

    if( integrate_exx ) {
      K = integrator.eval_exx(P, sn_link_settings);
      //matrix_type K_tmp = 0.5 * (K + K.transpose());
      //K = -K_tmp;
    } else { K = K_ref; }

#ifdef GAUXC_HAS_MPI
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    auto xc_int_end   = std::chrono::high_resolution_clock::now();
    double xc_int_dur = std::chrono::duration<double>( xc_int_end - xc_int_start ).count();

#ifdef GAUXC_HAS_MPI
    util::MPITimer mpi_lb_timings( MPI_COMM_WORLD, lb->get_timings() );
    util::MPITimer mpi_xc_timings( MPI_COMM_WORLD, integrator.get_timings() );
    util::MPITimer mpi_weight_timings( MPI_COMM_WORLD, mw.get_timings() );
#endif
    if( !world_rank ) {

      std::cout << std::scientific << std::setprecision(5) << std::endl;
      std::cout << "Load Balancer Timings" << std::endl;
      for( const auto& [name, dur] : lb->get_timings().all_timings() ) {
        #ifdef GAUXC_HAS_MPI
        const auto avg     = mpi_lb_timings.get_avg_duration(name).count();
        const auto min     = mpi_lb_timings.get_min_duration(name).count();
        const auto max     = mpi_lb_timings.get_max_duration(name).count();
        const auto std_dev = mpi_lb_timings.get_std_dev(name).count();
        #endif
        std::cout << "  " << std::setw(30) << name << ": " 
        #ifdef GAUXC_HAS_MPI
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
        #ifdef GAUXC_HAS_MPI
        const auto avg     = mpi_weight_timings.get_avg_duration(name).count();
        const auto min     = mpi_weight_timings.get_min_duration(name).count();
        const auto max     = mpi_weight_timings.get_max_duration(name).count();
        const auto std_dev = mpi_weight_timings.get_std_dev(name).count();
        #endif
        std::cout << "  " << std::setw(30) << name << ": " 
        #ifdef GAUXC_HAS_MPI
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
        #ifdef GAUXC_HAS_MPI
        const auto avg     = mpi_xc_timings.get_avg_duration(name).count();
        const auto min     = mpi_xc_timings.get_min_duration(name).count();
        const auto max     = mpi_xc_timings.get_max_duration(name).count();
        const auto std_dev = mpi_xc_timings.get_std_dev(name).count();
        #endif
        std::cout << "  " << std::setw(40) << name << ": " 
        #ifdef GAUXC_HAS_MPI
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
      if( not rks ) {
        std::cout << "| VXCz (ref)  |_F = " << VXCz_ref.norm() << std::endl;
        std::cout << "| VXCz (calc) |_F = " << VXCz.norm() << std::endl;
        std::cout << "RMS VXCz Diff     = " << (VXCz_ref - VXCz).norm() / basis.nbf()
                                           << std::endl;
      }
      if( gks ) {
        std::cout << "| VXCy (ref)  |_F = " << VXCy_ref.norm() << std::endl;
        std::cout << "| VXCy (calc) |_F = " << VXCy.norm() << std::endl;
        std::cout << "RMS VXCy Diff     = " << (VXCy_ref - VXCy).norm() / basis.nbf()
                                           << std::endl;
        std::cout << "| VXCx (ref)  |_F = " << VXCx_ref.norm() << std::endl;
        std::cout << "| VXCx (calc) |_F = " << VXCx.norm() << std::endl;
        std::cout << "RMS VXCx Diff     = " << (VXCx_ref - VXCx).norm() / basis.nbf()
                                           << std::endl;
      }
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

      if( not rks ) {
        dset = file.createDataSet<double>( "/DENSITY_Z", mat_space );
        dset.write_raw( Pz.data() );
      }
      if( gks ) {
        dset = file.createDataSet<double>( "/DENSITY_Y", mat_space );
        dset.write_raw( Py.data() );
        dset = file.createDataSet<double>( "/DENSITY_X", mat_space );
        dset.write_raw( Px.data() );
      }


      if( integrate_vxc ) {
        dset = file.createDataSet<double>( "/VXC", mat_space );
        dset.write_raw( VXC.data() );
        if( not rks ) {
          dset = file.createDataSet<double>( "/VXC_Z", mat_space );
          dset.write_raw( VXCz.data() );
        }
        if( gks ) {
          dset = file.createDataSet<double>( "/VXC_Y", mat_space );
          dset.write_raw( VXCy.data() );
          dset = file.createDataSet<double>( "/VXC_X", mat_space );
          dset.write_raw( VXCx.data() );
        }

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
#ifdef GAUXC_HAS_MPI
  MPI_Finalize();
#endif

}
