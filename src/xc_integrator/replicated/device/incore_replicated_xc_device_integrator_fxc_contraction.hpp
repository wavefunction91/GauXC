/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "incore_replicated_xc_device_integrator.hpp"
#include "vv10_nlc_device.hpp"
#include "device/xc_device_aos_data.hpp"
#include <gauxc/util/misc.hpp>
#include <gauxc/util/unused.hpp>

namespace GauXC::detail {

  template <typename ValueType>
  void IncoreReplicatedXCDeviceIntegrator<ValueType>::
    eval_fxc_contraction_( int64_t m, int64_t n, 
                          const value_type* P, int64_t ldp,
                          const value_type* tP, int64_t ldtp,
                          value_type* FXC, int64_t ldfxc,
                          const IntegratorSettingsXC& ks_settings ) {
    
    eval_fxc_contraction_( m, n, P, ldp, nullptr, 0, tP, ldtp, nullptr, 0,
                          FXC, ldfxc, nullptr, 0, ks_settings );
  }

    
  template <typename ValueType>
  void IncoreReplicatedXCDeviceIntegrator<ValueType>::
    eval_fxc_contraction_( int64_t m, int64_t n, 
                          const value_type* Ps, int64_t ldps,
                          const value_type* Pz, int64_t ldpz,
                          const value_type* tPs, int64_t ldtps,
                          const value_type* tPz, int64_t ldtpz,
                          value_type* FXCs, int64_t ldfxcs,
                          value_type* FXCz, int64_t ldfxcz,
                          const IntegratorSettingsXC& ks_settings ) {
    const bool is_uks = (Pz != nullptr);
    const bool is_rks = !is_uks;

    const auto& basis = this->load_balancer_->basis();

    // Check that P / FXC are sane
    const int64_t nbf = basis.nbf();
    if( m != n ) 
      GAUXC_GENERIC_EXCEPTION("P/FXC Must Be Square");
    if( m != nbf ) 
      GAUXC_GENERIC_EXCEPTION("P/FXC Must Have Same Dimension as Basis");
    if( ldps < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDPs");
    if( ldtps < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDTps");
    if( ldfxcs < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDFXCs");

    if( not is_rks ) {
      if( ldpz < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDPz");
      if( ldtpz < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDTpz");
      if( ldfxcz < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDFXCz");
    }

    // Get Tasks
    auto& tasks = this->load_balancer_->get_tasks();

    if( auto* nlc_settings = dynamic_cast<const IntegratorSettingsNLCInternal*>(&ks_settings) ) {
      if( not is_rks ) {
        GAUXC_GENERIC_EXCEPTION("NLC device FXC contraction currently supports only RKS");
      }
      if( this->reduction_driver_->takes_device_memory() ) {
        GAUXC_GENERIC_EXCEPTION("NLC device FXC contraction requires host reductions");
      }

      auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
      auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
      auto device_data_ptr = lwd->create_device_data(rt);

      GAUXC_MPI_CODE( MPI_Barrier(rt.comm());)

      value_type N_EL;
      this->timer_.time_op("XCIntegrator.LocalWork_FXC", [&](){
        nlc_fxc_contraction_local_work_( basis, Ps, ldps, tPs, ldtps, &N_EL,
          FXCs, ldfxcs, *nlc_settings, tasks.begin(), tasks.end(), *device_data_ptr );
      });

      GAUXC_MPI_CODE(
      this->timer_.time_op("XCIntegrator.ImbalanceWait_FXC",[&](){
        MPI_Barrier(this->load_balancer_->runtime().comm());
      });
      )

      this->timer_.time_op("XCIntegrator.Allreduce_FXC", [&](){
        this->reduction_driver_->allreduce_inplace( FXCs, nbf*nbf, ReductionOp::Sum );
        this->reduction_driver_->allreduce_inplace( &N_EL, 1, ReductionOp::Sum );
      });
      return;
    }

    // Allocate Device memory
    auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
    auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
    auto device_data_ptr = lwd->create_device_data(rt);

    GAUXC_MPI_CODE( MPI_Barrier(rt.comm());) 

    // Temporary electron count to judge integrator accuracy
    value_type N_EL;
  
    if( this->reduction_driver_->takes_device_memory() ) {

      // If we can do reductions on the device (e.g. NCCL)
      // Don't communicate data back to the host before reduction
      this->timer_.time_op("XCIntegrator.LocalWork_FXC", [&](){
        fxc_contraction_local_work_( basis, Ps, ldps, Pz, ldpz, tPs, ldtps, tPz, ldtpz,
          tasks.begin(), tasks.end(), *device_data_ptr);
      });

      GAUXC_MPI_CODE(
      this->timer_.time_op("XCIntegrator.ImbalanceWait_FXC",[&](){
        MPI_Barrier(this->load_balancer_->runtime().comm());
      });  
      )

      // Reduce results in device memory
      double* fxc_s_device = device_data_ptr->fxc_s_device_data();
      double* fxc_z_device;
      auto nel_device = device_data_ptr->nel_device_data();
      auto queue = device_data_ptr->queue();
      
      if( not is_rks ) {
        fxc_z_device = device_data_ptr->fxc_z_device_data();
        // UKS
        this->timer_.time_op("XCIntegrator.Allreduce_FXC", [&](){
          this->reduction_driver_->allreduce_inplace( fxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( fxc_z_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
        });
      } else {
        // RKS
        this->timer_.time_op("XCIntegrator.Allreduce_FXC", [&](){
          this->reduction_driver_->allreduce_inplace( fxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
        });
      }

      // Retrieve data to host
      this->timer_.time_op("XCIntegrator.DeviceToHostCopy_FXC",[&](){
        device_data_ptr->retrieve_fxc_contraction_integrands(&N_EL, FXCs, ldfxcs, FXCz, ldfxcz, nullptr, 0, nullptr, 0);
      });

    } else {

      // Compute local contributions to FXC and retrieve
      // data from device 
      this->timer_.time_op("XCIntegrator.LocalWork_FXC", [&](){
        fxc_contraction_local_work_( basis, Ps, ldps, Pz, ldpz, tPs, ldtps, tPz, ldtpz, &N_EL, 
                              FXCs, ldfxcs, FXCz, ldfxcz, tasks.begin(), tasks.end(), *device_data_ptr);
      });

      GAUXC_MPI_CODE(
      this->timer_.time_op("XCIntegrator.ImbalanceWait_FXC",[&](){
        MPI_Barrier(this->load_balancer_->runtime().comm());
      });  
      )

      // Reduce Results in host mem
      if( is_rks ) {
        this->timer_.time_op("XCIntegrator.Allreduce_FXC", [&](){
          this->reduction_driver_->allreduce_inplace( FXCs, nbf*nbf, ReductionOp::Sum );
        this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
        });
      } else {
        // UKS
        this->timer_.time_op("XCIntegrator.Allreduce_FXC", [&](){
          this->reduction_driver_->allreduce_inplace( FXCs, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( FXCz, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
        });
      }
    }
  }

  template <typename ValueType>
  void IncoreReplicatedXCDeviceIntegrator<ValueType>::
    nlc_fxc_contraction_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* tPs, int64_t ldtps,
                            value_type *N_EL,
                            value_type* FXC, int64_t ldfxc,
                            const IntegratorSettingsNLC& settings,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {

#ifdef GAUXC_HAS_CUDA
    auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
    auto* stack_data = dynamic_cast<XCDeviceStackData*>(&device_data);
    auto* aos_data = dynamic_cast<XCDeviceAoSData*>(&device_data);
    if( not stack_data or not aos_data ) GAUXC_BAD_LWD_DATA_CAST();

    const auto& func = *this->func_;
    if( not func.is_gga() ) {
      GAUXC_GENERIC_EXCEPTION("NLC device FXC contraction currently requires a GGA functional");
    }

    const auto& mol = this->load_balancer_->molecule();
    BasisSetMap basis_map(basis,mol);
    device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );

    auto task_comparator = []( const XCTask& a, const XCTask& b ) {
      return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
    };
    std::sort( task_begin, task_end, task_comparator );

    auto& lb_state = this->load_balancer_->state();
    if( not lb_state.modified_weights_are_stored ) {
      GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
    }

    integrator_term_tracker enabled_terms;
    enabled_terms.fxc_contraction = true;
    enabled_terms.ks_scheme = RKS;
    enabled_terms.xc_approx = integrator_xc_approx::GGA;

    const auto nbf = basis.nbf();
    const auto nshells = basis.nshells();
    device_data.reset_allocations();
    device_data.allocate_static_data_fxc_contraction( nbf, nshells, enabled_terms );
    device_data.send_static_data_density_basis( Ps, ldps, nullptr, 0, nullptr, 0, nullptr, 0, basis );
    device_data.send_static_data_trial_density( tPs, ldtps, nullptr, 0, nullptr, 0, nullptr, 0 );

    auto rt = detail::as_device_runtime(this->load_balancer_->runtime());
    auto* backend = rt.device_backend();

    std::vector<double> local_packed;
    auto task_it = task_begin;
    while( task_it != task_end ) {
      task_it = device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

      lwd->eval_collocation_gradient( &device_data );
      lwd->eval_xmat( 2.0, &device_data, false, DEN_S );
      lwd->eval_vvars_gga( &device_data, DEN_S );
      lwd->eval_uvars_gga( &device_data, RKS );
      lwd->eval_xmat_trial( 2.0, &device_data, false, DEN_S );
      lwd->eval_vvars_gga_trial( &device_data, DEN_S );

      const size_t batch_npts = stack_data->total_npts_task_batch;
      std::vector<double> points_x(batch_npts), points_y(batch_npts), points_z(batch_npts);
      std::vector<double> weights(batch_npts), rho(batch_npts), gamma(batch_npts);
      std::vector<double> grad_x(batch_npts), grad_y(batch_npts), grad_z(batch_npts);
      std::vector<double> trho(batch_npts), tgrad_x(batch_npts), tgrad_y(batch_npts), tgrad_z(batch_npts);

      auto base_stack = stack_data->base_stack;
      backend->copy_async( batch_npts, base_stack.points_x_device, points_x.data(), "NLC FXC points_x D2H" );
      backend->copy_async( batch_npts, base_stack.points_y_device, points_y.data(), "NLC FXC points_y D2H" );
      backend->copy_async( batch_npts, base_stack.points_z_device, points_z.data(), "NLC FXC points_z D2H" );
      backend->copy_async( batch_npts, base_stack.weights_device, weights.data(), "NLC FXC weights D2H" );
      backend->copy_async( batch_npts, base_stack.den_s_eval_device, rho.data(), "NLC FXC rho D2H" );
      backend->copy_async( batch_npts, base_stack.gamma_eval_device, gamma.data(), "NLC FXC gamma D2H" );
      backend->copy_async( batch_npts, base_stack.dden_sx_eval_device, grad_x.data(), "NLC FXC grad_x D2H" );
      backend->copy_async( batch_npts, base_stack.dden_sy_eval_device, grad_y.data(), "NLC FXC grad_y D2H" );
      backend->copy_async( batch_npts, base_stack.dden_sz_eval_device, grad_z.data(), "NLC FXC grad_z D2H" );
      backend->copy_async( batch_npts, base_stack.tden_s_eval_device, trho.data(), "NLC FXC trho D2H" );
      backend->copy_async( batch_npts, base_stack.tdden_sx_eval_device, tgrad_x.data(), "NLC FXC tgrad_x D2H" );
      backend->copy_async( batch_npts, base_stack.tdden_sy_eval_device, tgrad_y.data(), "NLC FXC tgrad_y D2H" );
      backend->copy_async( batch_npts, base_stack.tdden_sz_eval_device, tgrad_z.data(), "NLC FXC tgrad_z D2H" );
      backend->master_queue_synchronize();

      local_packed.reserve(local_packed.size() + 13 * batch_npts);
      for( size_t i = 0; i < batch_npts; ++i ) {
        local_packed.push_back(points_x[i]);
        local_packed.push_back(points_y[i]);
        local_packed.push_back(points_z[i]);
        local_packed.push_back(weights[i]);
        local_packed.push_back(rho[i]);
        local_packed.push_back(gamma[i]);
        local_packed.push_back(grad_x[i]);
        local_packed.push_back(grad_y[i]);
        local_packed.push_back(grad_z[i]);
        local_packed.push_back(trho[i]);
        local_packed.push_back(tgrad_x[i]);
        local_packed.push_back(tgrad_y[i]);
        local_packed.push_back(tgrad_z[i]);
      }
    }

    size_t local_point_offset = 0;
    const auto global_packed = vv10::allgather_packed_grid(
      *this->reduction_driver_, local_packed, 13, local_point_offset );
    if( global_packed.size() % 13 != 0 ) {
      GAUXC_GENERIC_EXCEPTION("Invalid VV10 device FXC packed grid size after allgather");
    }

    const size_t global_npts = global_packed.size() / 13;
    std::vector<double> coords(3 * global_npts), weights(global_npts), rho(global_npts), gamma(global_npts);
    std::vector<double> grad_x(global_npts), grad_y(global_npts), grad_z(global_npts);
    std::vector<double> trho(global_npts), tgrad_x(global_npts), tgrad_y(global_npts), tgrad_z(global_npts);
    std::vector<double> fxc_A(global_npts, 0.0), fxc_B(3 * global_npts, 0.0);
    for( size_t i = 0; i < global_npts; ++i ) {
      coords[3*i+0] = global_packed[13*i+0];
      coords[3*i+1] = global_packed[13*i+1];
      coords[3*i+2] = global_packed[13*i+2];
      weights[i] = global_packed[13*i+3];
      rho[i] = global_packed[13*i+4];
      gamma[i] = global_packed[13*i+5];
      grad_x[i] = global_packed[13*i+6];
      grad_y[i] = global_packed[13*i+7];
      grad_z[i] = global_packed[13*i+8];
      trho[i] = global_packed[13*i+9];
      tgrad_x[i] = global_packed[13*i+10];
      tgrad_y[i] = global_packed[13*i+11];
      tgrad_z[i] = global_packed[13*i+12];
    }

    vv10::eval_fxc_rks_cuda( 0, settings,
      vv10::ResponseGridView{ global_npts, coords.data(), weights.data(), rho.data(), gamma.data(),
        grad_x.data(), grad_y.data(), grad_z.data() },
      vv10::TrialView{ trho.data(), tgrad_x.data(), tgrad_y.data(), tgrad_z.data() },
      vv10::ResponseCorrectionsView{ fxc_A.data(), fxc_B.data() } );

    device_data.zero_fxc_contraction_integrands();
    task_it = task_begin;
    size_t local_batch_offset = 0;
    while( task_it != task_end ) {
      task_it = device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

      lwd->eval_collocation_gradient( &device_data );
      lwd->eval_xmat( 2.0, &device_data, false, DEN_S );
      lwd->eval_vvars_gga( &device_data, DEN_S );
      lwd->eval_uvars_gga( &device_data, RKS );
      lwd->inc_nel( &device_data );

      const size_t batch_npts = stack_data->total_npts_task_batch;
      const size_t global_offset = local_point_offset + local_batch_offset;
      auto base_stack = stack_data->base_stack;
      std::vector<double> batch_A(batch_npts), batch_Bx(batch_npts), batch_By(batch_npts), batch_Bz(batch_npts);

      for( size_t i = 0; i < batch_npts; ++i ) {
        const auto global_i = global_offset + i;
        batch_A[i]  = weights[global_i] * fxc_A[global_i];
        batch_Bx[i] = weights[global_i] * fxc_B[3*global_i+0];
        batch_By[i] = weights[global_i] * fxc_B[3*global_i+1];
        batch_Bz[i] = weights[global_i] * fxc_B[3*global_i+2];
      }

      backend->copy_async( batch_npts, batch_A.data(), base_stack.FXC_A_s_eval_device, "NLC FXC A H2D" );
      backend->copy_async( batch_npts, batch_Bx.data(), base_stack.FXC_Bx_s_eval_device, "NLC FXC Bx H2D" );
      backend->copy_async( batch_npts, batch_By.data(), base_stack.FXC_By_s_eval_device, "NLC FXC By H2D" );
      backend->copy_async( batch_npts, batch_Bz.data(), base_stack.FXC_Bz_s_eval_device, "NLC FXC Bz H2D" );
      backend->master_queue_synchronize();

      lwd->eval_zmat_gga_fxc( &device_data, DEN_S );
      lwd->inc_fxc( &device_data, DEN_S, false );
      local_batch_offset += batch_npts;
    }

    lwd->symmetrize_fxc( &device_data, DEN_S );
    backend->master_queue_synchronize();
    device_data.retrieve_fxc_contraction_integrands( N_EL, FXC, ldfxc, nullptr, 0, nullptr, 0, nullptr, 0 );
#else
    GAUXC_GENERIC_EXCEPTION("NLC device FXC contraction requires CUDA");
#endif
  }

  template <typename ValueType>
  void IncoreReplicatedXCDeviceIntegrator<ValueType>::
    fxc_contraction_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* tPs, int64_t ldtps,
                            const value_type* tPz, int64_t ldtpz,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data) {
    const bool is_uks = (Pz != nullptr);
    const bool is_rks = !is_uks;
    if (not is_rks and not is_uks) {
      GAUXC_GENERIC_EXCEPTION("MUST BE EITHER RKS OR UKS!");
    }
    

    // Cast LWD to LocalDeviceWorkDriver
    auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

    // Setup Aliases
    const auto& func  = *this->func_;
    const auto& mol   = this->load_balancer_->molecule();

    // Get basis map
    BasisSetMap basis_map(basis,mol);

    // Populate submat maps
    device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );


    // Sort tasks 
    auto task_comparator = []( const XCTask& a, const XCTask& b ) {
      return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
    };
    std::sort( task_begin, task_end, task_comparator );


    // Check that Partition Weights have been calculated
    auto& lb_state = this->load_balancer_->state();
    if( not lb_state.modified_weights_are_stored ) {
      GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified"); 
    }
    

    integrator_term_tracker enabled_terms;
    enabled_terms.fxc_contraction = true;

    if (is_rks) enabled_terms.ks_scheme = RKS;
    else if (is_uks) enabled_terms.ks_scheme = UKS;

    if( func.is_lda() )      
      enabled_terms.xc_approx = integrator_xc_approx::LDA; 
    else if( func.is_gga() ) 
      enabled_terms.xc_approx = integrator_xc_approx::GGA; 
    else if( func.needs_laplacian() )                    
      GAUXC_GENERIC_EXCEPTION("FXC contraction does not support MGGA with Laplacian");
    else
      enabled_terms.xc_approx = integrator_xc_approx::MGGA_TAU;
    
    // Do XC integration in task batches
    const auto nbf     = basis.nbf();
    const auto nshells = basis.nshells();
    device_data.reset_allocations();
    device_data.allocate_static_data_fxc_contraction( nbf, nshells, enabled_terms);
    
    device_data.send_static_data_density_basis( Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, basis );
    device_data.send_static_data_trial_density( tPs, ldtps, tPz, ldtpz, nullptr, 0, nullptr, 0 );


    // Zero integrands
    device_data.zero_fxc_contraction_integrands();


    auto task_it = task_begin;
    while( task_it != task_end ) {

      // Determine next task batch, send relevant data to device (FXC only)
      task_it = 
        device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

      /*** Process the batches ***/
      
      const bool need_lapl = func.needs_laplacian();
      // Evaluate collocation
      if( func.is_mgga() ) {
        if(need_lapl) lwd->eval_collocation_laplacian( &device_data );
        else          lwd->eval_collocation_gradient( &device_data );
      }
      else if( func.is_gga() ) lwd->eval_collocation_gradient( &device_data );
      else                     lwd->eval_collocation( &device_data );
        
      const double xmat_fac = is_rks ? 2.0 : 1.0;
      const bool need_xmat_grad = func.is_mgga();

      // Evaluate X matrix and V vars
      auto do_xmat_vvar = [&](density_id den_id) {
        lwd->eval_xmat( xmat_fac, &device_data, need_xmat_grad, den_id );
        if(func.is_lda())      lwd->eval_vvars_lda( &device_data, den_id );
        else if(func.is_gga()) lwd->eval_vvars_gga( &device_data, den_id ); 
        else                   lwd->eval_vvars_mgga( &device_data, den_id, need_lapl );
      };

      do_xmat_vvar(DEN_S);
      if (not is_rks) {
        do_xmat_vvar(DEN_Z);
      }

      // Evaluate U variables
      if( func.is_mgga() )      lwd->eval_uvars_mgga( &device_data, enabled_terms.ks_scheme, need_lapl );
      else if( func.is_gga() )  lwd->eval_uvars_gga ( &device_data, enabled_terms.ks_scheme );
      else                      lwd->eval_uvars_lda ( &device_data, enabled_terms.ks_scheme );

      // Evaluate XC functional
      if( func.is_mgga() )     lwd->eval_kern_vxc_fxc_mgga( func, &device_data );
      else if( func.is_gga() ) lwd->eval_kern_vxc_fxc_gga ( func, &device_data );
      else                     lwd->eval_kern_vxc_fxc_lda ( func, &device_data );      

      // Do scalar N_EL integrations
      lwd->inc_nel( &device_data );

      
      // Evaluate X matrix and V vars from trial density
      auto do_xmat_vvar_trial = [&](density_id den_id) {
        lwd->eval_xmat_trial( xmat_fac, &device_data, need_xmat_grad, den_id );
        if(func.is_lda())      lwd->eval_vvars_lda_trial( &device_data, den_id );
        else if(func.is_gga()) lwd->eval_vvars_gga_trial( &device_data, den_id ); 
        else                   lwd->eval_vvars_mgga_trial( &device_data, den_id, need_lapl );
      };

      do_xmat_vvar_trial(DEN_S);
      if (not is_rks) {
        do_xmat_vvar_trial(DEN_Z);
      }

      // Evaluate tmat (it contains the trial u variable evaluation inside)
      if( func.is_mgga() )      lwd->eval_tmat_mgga( &device_data, enabled_terms.ks_scheme, need_lapl );
      else if( func.is_gga() )  lwd->eval_tmat_gga ( &device_data, enabled_terms.ks_scheme );
      else                      lwd->eval_tmat_lda ( &device_data, enabled_terms.ks_scheme );

      auto do_zmat_fxc = [&](density_id den_id) {
        if( func.is_mgga() ) {
          lwd->eval_zmat_mgga_fxc( &device_data, need_lapl, den_id);
          lwd->eval_mmat_mgga_fxc( &device_data, need_lapl, den_id);
        }
        else if( func.is_gga() ) 
          lwd->eval_zmat_gga_fxc( &device_data, den_id );
        else 
          lwd->eval_zmat_lda_fxc( &device_data, den_id );
        lwd->inc_fxc( &device_data, den_id, func.is_mgga() );
      };

      do_zmat_fxc(DEN_S);
      if(not is_rks) {
        do_zmat_fxc(DEN_Z);
      } 

    } // Loop over batches of batches 

    // Symmetrize FXC in device memory
    lwd->symmetrize_fxc( &device_data, DEN_S );
    if (not is_rks) {
      lwd->symmetrize_fxc( &device_data, DEN_Z );
      }
  }

  template <typename ValueType>
  void IncoreReplicatedXCDeviceIntegrator<ValueType>::
    fxc_contraction_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* tPs, int64_t ldtps,
                            const value_type* tPz, int64_t ldtpz,
                            value_type *N_EL,
                            value_type* FXCs, int64_t ldfxcs,
                            value_type* FXCz, int64_t ldfxcz,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {
    
    // Get integrate and keep data on device
    fxc_contraction_local_work_( basis, Ps, ldps, Pz, ldpz, tPs, ldtps, tPz, ldtpz, 
                              task_begin, task_end, device_data);
    auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
    rt.device_backend()->master_queue_synchronize();

    // Receive FXC terms from host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_FXC",[&](){
      device_data.retrieve_fxc_contraction_integrands( N_EL, FXCs, ldfxcs, FXCz, ldfxcz, nullptr, 0, nullptr, 0 ); 
    });
  }
}
