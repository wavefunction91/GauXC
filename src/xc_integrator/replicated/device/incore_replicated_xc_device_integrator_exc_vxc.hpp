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
#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include "device/xc_device_stack_data.hpp"
#include "vv10_nlc_device.hpp"
#include <fstream>
#include <gauxc/exceptions.hpp>
#include <gauxc/util/unused.hpp>
#include <numeric>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC, const IntegratorSettingsXC& settings ) {
  eval_exc_vxc_( m, n, P, ldp, nullptr, 0, nullptr, 0, nullptr, 0, 
                      VXC, ldvxc, nullptr, 0, nullptr, 0, nullptr, 0, EXC, settings );
}


template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC, const IntegratorSettingsXC& settings ) { 
  eval_exc_vxc_( m, n, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, 
                VXCs, ldvxcs, VXCz, ldvxcz, nullptr, 0, nullptr, 0, EXC, settings );
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      const value_type* Py,
                      int64_t ldpy,
                      const value_type* Px,
                      int64_t ldpx,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* VXCy, int64_t ldvxcy,
                      value_type* VXCx, int64_t ldvxcx,
                      value_type* EXC, const IntegratorSettingsXC& settings ) {
  const bool is_gks = (Pz != nullptr) and (Py != nullptr) and (Px != nullptr);
  const bool is_uks = (Pz != nullptr) and (Py == nullptr) and (Px == nullptr);
  const bool is_rks = (Ps != nullptr) and (not is_uks and not is_gks);
  const auto* nlc_settings = dynamic_cast<const IntegratorSettingsNLCInternal*>(&settings);

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldps < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPs");
  if( ldvxcs < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCs");

  if( not is_rks ) {
    if( ldpz < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDPz");
    if( ldvxcz < nbf )
      GAUXC_GENERIC_EXCEPTION("Invalid LDVXCz");
    if( is_gks ) {
      if( ldpy < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDPy");
      if( ldvxcy < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDVXCy");
      if( ldpx < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDPx");
      if( ldvxcx < nbf )
        GAUXC_GENERIC_EXCEPTION("Invalid LDVXCx");
    }
  }

  if( nlc_settings ) {
    if( is_gks ) {
      GAUXC_GENERIC_EXCEPTION("NLC device EXC/VXC is not implemented for GKS");
    }
    if( this->reduction_driver_->takes_device_memory() ) {
      GAUXC_GENERIC_EXCEPTION("NLC device EXC/VXC currently requires host reductions");
    }

    auto& tasks = this->load_balancer_->get_tasks();
    auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
    auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
    auto device_data_ptr = lwd->create_device_data(rt);
    value_type N_EL;
    if( VXCz ) std::fill_n( VXCz, nbf*nbf, 0.0 );

    this->timer_.time_op("XCIntegrator.LocalWork_NLC_EXC_VXC", [&](){
      nlc_exc_vxc_local_work_( basis, Ps, ldps, is_rks ? 2.0 : 1.0,
        VXCs, ldvxcs, EXC, &N_EL,
        *nlc_settings, tasks.begin(), tasks.end(), *device_data_ptr );
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_NLC_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });
    )

    this->timer_.time_op("XCIntegrator.Allreduce_NLC_EXC_VXC", [&](){
      if( VXCs ) this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
      if( VXCz ) this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( EXC, 1, ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( &N_EL, 1, ReductionOp::Sum );
    });
    return;
  }


  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

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
    this->timer_.time_op("XCIntegrator.LocalWork_EXC_VXC", [&](){
      exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, Py, ldpy, Px, ldpx, tasks.begin(), tasks.end(), 
        *device_data_ptr, true);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce results in device memory
    double* vxc_s_device = device_data_ptr->vxc_s_device_data();
    double* vxc_z_device;
    double* vxc_y_device;
    double* vxc_x_device;
    auto exc_device = device_data_ptr->exc_device_data();
    auto nel_device = device_data_ptr->nel_device_data();
    auto queue      = device_data_ptr->queue();
    
    if( not is_rks ) {
      vxc_z_device = device_data_ptr->vxc_z_device_data();
      if( is_gks ) {
        // GKS
        vxc_y_device = device_data_ptr->vxc_y_device_data();
        vxc_x_device = device_data_ptr->vxc_x_device_data();
        this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
          this->reduction_driver_->allreduce_inplace( vxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( vxc_z_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( vxc_y_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( vxc_x_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
        });
      } else {
        // UKS
        this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
          this->reduction_driver_->allreduce_inplace( vxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( vxc_z_device, nbf*nbf, ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
          this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
        });

      }
    } else {
      // RKS
      this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
        this->reduction_driver_->allreduce_inplace( vxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
        this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
        this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
      });
    }


    // Retrieve data to host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
      device_data_ptr->retrieve_exc_vxc_integrands( EXC, &N_EL, VXCs, ldvxcs, VXCz, ldvxcz,
                                                                VXCy, ldvxcy, VXCx, ldvxcx );
    });


  } else {

    // Compute local contributions to EXC/VXC and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork_EXC_VXC", [&](){
      exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, Py, ldpy, Px, ldpx,
                                VXCs, ldvxcs, VXCz, ldvxcz, VXCy, ldvxcy, VXCx, ldvxcx, EXC, 
                              &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce Results in host mem
    if( is_rks ) {
      this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
        this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
        this->reduction_driver_->allreduce_inplace( EXC, 1,       ReductionOp::Sum );
        this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
      });
    } else {
      if( is_gks ) {
        this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
          this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( VXCy, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( VXCx, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( EXC, 1,       ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
        });
      } else {
        // UKS
        this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
          this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( EXC, 1,       ReductionOp::Sum );
          this->reduction_driver_->allreduce_inplace( &N_EL, 1,       ReductionOp::Sum );
        });

      }
    }
  }
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  nlc_exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            double xmat_fac,
                            value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                            const IntegratorSettingsNLC& settings,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto* stack_data = dynamic_cast<XCDeviceStackData*>(&device_data);
  if( not stack_data ) GAUXC_BAD_LWD_DATA_CAST();

  const auto& func = *this->func_;
  const auto& mol = this->load_balancer_->molecule();
  if( not func.is_gga() ) {
    GAUXC_GENERIC_EXCEPTION("NLC device EXC/VXC currently requires a GGA functional");
  }

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
  enabled_terms.exc_vxc = true;
  enabled_terms.ks_scheme = RKS;
  enabled_terms.xc_approx = integrator_xc_approx::GGA;

  const auto nbf = basis.nbf();
  const auto nshells = basis.nshells();
  device_data.reset_allocations();
  device_data.allocate_static_data_exc_vxc( nbf, nshells, enabled_terms, true );
  device_data.send_static_data_density_basis( Ps, ldps, nullptr, 0, nullptr, 0, nullptr, 0, basis );

  size_t local_npts = 0;
  const size_t ntasks = std::distance(task_begin, task_end);
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    local_npts += (task_begin + iT)->points.size();
  }

  std::vector<double> local_packed;
  local_packed.reserve(6 * local_npts);

  auto rt = detail::as_device_runtime(this->load_balancer_->runtime());
  auto* backend = rt.device_backend();

  auto task_it = task_begin;
  while( task_it != task_end ) {
    const auto batch_begin = task_it;
    task_it = device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    lwd->eval_collocation_gradient( &device_data );
    lwd->eval_xmat( xmat_fac, &device_data, false, DEN_S );
    lwd->eval_vvars_gga( &device_data, DEN_S );
    lwd->eval_uvars_gga( &device_data, RKS );

    const size_t batch_npts = stack_data->total_npts_task_batch;
    std::vector<double> rho(batch_npts), gamma(batch_npts);
    auto base_stack = stack_data->base_stack;
    backend->copy_async( batch_npts, base_stack.den_s_eval_device, rho.data(), "NLC rho D2H" );
    backend->copy_async( batch_npts, base_stack.gamma_eval_device, gamma.data(), "NLC gamma D2H" );
    backend->master_queue_synchronize();

    size_t batch_offset = 0;
    for( auto batch_task = batch_begin; batch_task != task_it; ++batch_task ) {
      const auto npts = batch_task->points.size();
      const auto* points = batch_task->points.data()->data();
      const auto* weights = batch_task->weights.data();
      for( size_t i = 0; i < npts; ++i, ++batch_offset ) {
        local_packed.push_back(points[3*i + 0]);
        local_packed.push_back(points[3*i + 1]);
        local_packed.push_back(points[3*i + 2]);
        local_packed.push_back(weights[i]);
        local_packed.push_back(rho[batch_offset]);
        local_packed.push_back(gamma[batch_offset]);
      }
    }
    if( batch_offset != batch_npts ) {
      GAUXC_GENERIC_EXCEPTION("Invalid NLC device EXC/VXC host/device batch point count");
    }
  }

  size_t local_point_offset = 0;
  const auto global_packed = vv10::allgather_packed_grid(
    *this->reduction_driver_, local_packed, 6, local_point_offset );
  if( global_packed.size() % 6 != 0 ) {
    GAUXC_GENERIC_EXCEPTION("Invalid VV10 device EXC/VXC packed grid size after allgather");
  }

  const size_t global_npts = global_packed.size() / 6;
  std::vector<double> coords(3 * global_npts), weights(global_npts), rho(global_npts), gamma(global_npts);
  std::vector<double> eps(global_npts, 0.0), vrho(global_npts, 0.0), vgamma(global_npts, 0.0);
  for( size_t i = 0; i < global_npts; ++i ) {
    coords[3*i+0] = global_packed[6*i+0];
    coords[3*i+1] = global_packed[6*i+1];
    coords[3*i+2] = global_packed[6*i+2];
    weights[i] = global_packed[6*i+3];
    rho[i] = global_packed[6*i+4];
    gamma[i] = global_packed[6*i+5];
  }

  vv10::eval_exc_vxc_cuda( 0, settings,
    vv10::GridView{ global_npts, coords.data(), weights.data(), rho.data(), gamma.data() },
    vv10::CorrectionsView{ eps.data(), vrho.data(), vgamma.data() } );

  device_data.zero_exc_vxc_integrands(enabled_terms);

  task_it = task_begin;
  size_t local_batch_offset = 0;
  while( task_it != task_end ) {
    task_it = device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    lwd->eval_collocation_gradient( &device_data );
    lwd->eval_xmat( xmat_fac, &device_data, false, DEN_S );
    lwd->eval_vvars_gga( &device_data, DEN_S );
    lwd->eval_uvars_gga( &device_data, RKS );

    const size_t batch_npts = stack_data->total_npts_task_batch;
  const size_t global_offset = local_point_offset + local_batch_offset;

    auto base_stack = stack_data->base_stack;
    std::vector<double> batch_eps(batch_npts), batch_vrho(batch_npts), batch_vgamma(batch_npts);
    for( size_t i = 0; i < batch_npts; ++i ) {
      const auto global_i = global_offset + i;
      batch_eps[i] = eps[global_i] * weights[global_i];
      batch_vrho[i] = vrho[global_i] * weights[global_i];
      batch_vgamma[i] = vgamma[global_i] * weights[global_i];
    }
    backend->copy_async( batch_npts, batch_eps.data(), base_stack.eps_eval_device, "NLC eps H2D" );
    backend->copy_async( batch_npts, batch_vrho.data(), base_stack.vrho_eval_device, "NLC vrho H2D" );
    backend->copy_async( batch_npts, batch_vgamma.data(), base_stack.vgamma_eval_device, "NLC vgamma H2D" );
    backend->master_queue_synchronize();

    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );
    if( VXC ) {
      lwd->eval_zmat_gga_vxc( &device_data, RKS, DEN_S );
      lwd->inc_vxc( &device_data, DEN_S, false );
    }
    local_batch_offset += batch_npts;
  }

  if( VXC ) lwd->symmetrize_vxc( &device_data, DEN_S );
  backend->master_queue_synchronize();
  device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXC, ldvxc, nullptr, 0, nullptr, 0, nullptr, 0 );
}


template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data, bool do_vxc ) {
  const bool is_gks = (Pz != nullptr) and (Py != nullptr) and (Px != nullptr);
  const bool is_uks = (Pz != nullptr) and (Py == nullptr) and (Px == nullptr);
  const bool is_rks = (Ps != nullptr) and (not is_uks and not is_gks);
  if (not is_rks and not is_uks and not is_gks) {
    GAUXC_GENERIC_EXCEPTION("MUST BE EITHER RKS, UKS, or GKS!");
  }
  

  // Cast LWD to LocalDeviceWorkDriver
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& mol   = this->load_balancer_->molecule();

  if( func.is_mgga() and is_gks ) GAUXC_GENERIC_EXCEPTION("GKS mGGAs NYI!");

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
  enabled_terms.exc_vxc = true;

  if (is_rks) enabled_terms.ks_scheme = RKS;
  else if (is_uks) enabled_terms.ks_scheme = UKS;
  else if (is_gks) enabled_terms.ks_scheme = GKS;

  if( func.is_lda() )      
    enabled_terms.xc_approx = integrator_xc_approx::LDA; 
  else if( func.is_gga() ) 
    enabled_terms.xc_approx = integrator_xc_approx::GGA; 
  else if( func.needs_laplacian() )                    
    enabled_terms.xc_approx = integrator_xc_approx::MGGA_LAPL;
  else
    enabled_terms.xc_approx = integrator_xc_approx::MGGA_TAU;
  
  // Do XC integration in task batches
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  device_data.reset_allocations();
  device_data.allocate_static_data_exc_vxc( nbf, nshells, enabled_terms, do_vxc );
  
  device_data.send_static_data_density_basis( Ps, ldps, Pz, ldpz, Px, ldpx, Py, ldpy, basis );



  // Zero integrands
  device_data.zero_exc_vxc_integrands(enabled_terms);


  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXC VXC only)
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
      if (not is_uks) {
        do_xmat_vvar(DEN_Y);
        do_xmat_vvar(DEN_X);
      }
    }


    // Evaluate U variables
    if( func.is_mgga() )      lwd->eval_uvars_mgga( &device_data, enabled_terms.ks_scheme, need_lapl );
    else if( func.is_gga() )  lwd->eval_uvars_gga ( &device_data, enabled_terms.ks_scheme );
    else                      lwd->eval_uvars_lda ( &device_data, enabled_terms.ks_scheme );

    // Evaluate XC functional
    if( func.is_mgga() )     lwd->eval_kern_exc_vxc_mgga( func, &device_data );
    else if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga ( func, &device_data );
    else                     lwd->eval_kern_exc_vxc_lda ( func, &device_data );
    

    // Do scalar EXC/N_EL integrations
    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );
    if( not do_vxc) continue;

   auto do_zmat_vxc = [&](density_id den_id) {
     if( func.is_mgga() ) {
       lwd->eval_zmat_mgga_vxc( &device_data, enabled_terms.ks_scheme, need_lapl, den_id);
       lwd->eval_mmat_mgga_vxc( &device_data, enabled_terms.ks_scheme, need_lapl, den_id);
     }
     else if( func.is_gga() ) 
       lwd->eval_zmat_gga_vxc( &device_data, enabled_terms.ks_scheme, den_id );
     else 
       lwd->eval_zmat_lda_vxc( &device_data, enabled_terms.ks_scheme, den_id );
     lwd->inc_vxc( &device_data, den_id, func.is_mgga() );
  };

  do_zmat_vxc(DEN_S);
  if(not is_rks) {
    do_zmat_vxc(DEN_Z);
    if(not is_uks) {
      do_zmat_vxc(DEN_Y);
      do_zmat_vxc(DEN_X);
    }
  } 

  } // Loop over batches of batches 

  // Symmetrize VXC in device memory
  if( do_vxc ) {
    lwd->symmetrize_vxc( &device_data, DEN_S );
    if (not is_rks) {
      lwd->symmetrize_vxc( &device_data, DEN_Z );
      if (not is_uks) {
        lwd->symmetrize_vxc( &device_data, DEN_Y );
        lwd->symmetrize_vxc( &device_data, DEN_X );
      }
    }
  }
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,   
                            value_type* VXCs, int64_t ldvxcs,
                            value_type* VXCz, int64_t ldvxcz,
                            value_type* VXCy, int64_t ldvxcy,
                            value_type* VXCx, int64_t ldvxcx, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {
  
  // Get integrate and keep data on device
  const bool do_vxc = VXCs;
  exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, Py, ldpy, Px, ldpx, task_begin, task_end, device_data, do_vxc );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  rt.device_backend()->master_queue_synchronize();

  // Receive XC terms from host
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
    device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXCs, ldvxcs, VXCz, ldvxcz, VXCy, ldvxcy, VXCx, ldvxcx ); 
  });

}


}
}

