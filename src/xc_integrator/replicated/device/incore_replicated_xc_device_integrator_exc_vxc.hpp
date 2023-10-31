/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/exceptions.hpp>
#include <gauxc/util/unused.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC, const IntegratorSettingsXC& settings ) {


  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxc < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");


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
    // Don't communicate data back to the hot before reduction
    this->timer_.time_op("XCIntegrator.LocalWork_EXC_VXC", [&](){
      exc_vxc_local_work_( basis, P, ldp, tasks.begin(), tasks.end(), 
        *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce results in device memory
    auto vxc_device = device_data_ptr->vxc_device_data();
    auto exc_device = device_data_ptr->exc_device_data();
    auto nel_device = device_data_ptr->nel_device_data();
    auto queue      = device_data_ptr->queue();
    this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
      this->reduction_driver_->allreduce_inplace( vxc_device, nbf*nbf, ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
    });

    // Retrieve data to host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
      device_data_ptr->retrieve_exc_vxc_integrands( EXC, &N_EL, VXC, ldvxc );
    });


  } else {

    // Compute local contributions to EXC/VXC and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork_EXC_VXC", [&](){
      exc_vxc_local_work_( basis, P, ldp, VXC, ldvxc, EXC, 
        &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce Results in host mem
    this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
      this->reduction_driver_->allreduce_inplace( VXC, nbf*nbf, ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );
    });

  }
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
  GauXC::util::unused(m,n,Ps,ldps,Pz,ldpz,VXCs,ldvxcs,VXCz,ldvxcz,EXC,settings);
  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");
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
  GauXC::util::unused(m,n,Ps,ldps,Pz,ldpz,Py,ldpy,Px,ldpx,VXCs,ldvxcs,VXCz,ldvxcz,VXCy,ldvxcy,VXCx,ldvxcx,EXC,settings);
  GAUXC_GENERIC_EXCEPTION("GKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  neo_eval_exc_vxc_( int64_t m1, int64_t n1, int64_t m2, int64_t n2, 
                    const value_type* P1s, int64_t ldp1s,
                    const value_type* P2s, int64_t ldp2s,
                    const value_type* P2z, int64_t ldp2z,
                    value_type* VXC1s, int64_t ldvxc1s,
                    value_type* VXC2s, int64_t ldvxc2s,
                    value_type* VXC2z, int64_t ldvxc2z,
                    value_type* EXC1,  value_type* EXC2, const IntegratorSettingsXC& settings ) {
  GauXC::util::unused(m1,n1,m2,n2,P1s,ldp1s,P2s,ldp2s,P2z,ldp2z,VXC1s,ldvxc1s,VXC2s,ldvxc2s,VXC2z,ldvxc2z,EXC1,EXC2,settings);
  GAUXC_GENERIC_EXCEPTION("NEO-RKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_neo_exc_vxc_( int64_t m1, int64_t n1, int64_t m2, int64_t n2, 
                    const value_type* P1s, int64_t ldp1s,
                    const value_type* P1z, int64_t ldp1z,
                    const value_type* P2s, int64_t ldp2s,
                    const value_type* P2z, int64_t ldp2z,
                    value_type* VXC1s, int64_t ldvxc1s,
                    value_type* VXC1z, int64_t ldvxc1z,
                    value_type* VXC2s, int64_t ldvxc2s,
                    value_type* VXC2z, int64_t ldvxc2z,
                    value_type* EXC1,  value_type* EXC2, const IntegratorSettingsXC& settings ) {

  GauXC::util::unused(m1,n1,m2,n2,P1s,ldp1s,P1z,ldp1z,P2s,ldp2s,P2z,ldp2z,VXC1s,ldvxc1s,VXC1z,ldvxc1z,VXC2s,ldvxc2s,VXC2z,ldvxc2z,EXC1,EXC2,settings);
  GAUXC_GENERIC_EXCEPTION("NEO-UKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data ) {


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
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified"); 
  }

#if 0
  this->timer_.time_op("XCIntegrator.ScreenWeights",[&](){

  constexpr double weight_thresh = std::numeric_limits<double>::epsilon();
  for( auto it = task_begin; it != task_end; ++it ) {
    it->max_weight = *std::max_element( it->weights.begin(), it->weights.end() );
  }

  size_t old_ntasks = std::distance( task_begin, task_end );
  task_end = std::stable_partition(task_begin, task_end,
    [&](const auto& a){ return a.max_weight > weight_thresh; } );

  size_t new_ntasks = std::distance( task_begin, task_end );
  std::cout << old_ntasks << ", " << new_ntasks << std::endl;

  });
#endif

  // Do XC integration in task batches
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  device_data.reset_allocations();
  device_data.allocate_static_data_exc_vxc( nbf, nshells );
  device_data.send_static_data_density_basis( P, ldp, basis );

  // Zero integrands
  device_data.zero_exc_vxc_integrands();

  // Processes batches in groups that saturadate available device memory
  integrator_term_tracker enabled_terms;
  enabled_terms.exc_vxc = true;
  if( func.is_lda() )      
    enabled_terms.xc_approx = integrator_xc_approx::LDA; 
  else if( func.is_gga() ) 
    enabled_terms.xc_approx = integrator_xc_approx::GGA; 
  else if( func.needs_laplacian() )                    
    enabled_terms.xc_approx = integrator_xc_approx::MGGA_LAPL;
  else
    enabled_terms.xc_approx = integrator_xc_approx::MGGA_TAU;

#if 0
  std::vector<int32_t> full_shell_list(nshells);
  std::iota(full_shell_list.begin(),full_shell_list.end(),0);
  for( auto it = task_begin; it != task_end; ++it ) {
    it->cou_screening.shell_list = full_shell_list;
    it->cou_screening.nbe        = nbf;
  }
#endif

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

    // Evaluate X matrix
    const bool need_xmat_grad = func.is_mgga();
    lwd->eval_xmat( 2.0, &device_data, need_xmat_grad );

    // Evaluate U/V variables
    if( func.is_mgga() )     lwd->eval_uvvar_mgga_rks( &device_data, need_lapl );
    else if( func.is_gga() ) lwd->eval_uvvar_gga_rks( &device_data );
    else                     lwd->eval_uvvar_lda_rks( &device_data );

    // Evaluate XC functional
    if( func.is_mgga() )     lwd->eval_kern_exc_vxc_mgga( func, &device_data );
    else if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                     lwd->eval_kern_exc_vxc_lda( func, &device_data );

    // Do scalar EXC/N_EL integrations
    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );

    // Evaluate Z (+ M) matrix
    if( func.is_mgga() ) {
      lwd->eval_zmat_mgga_vxc_rks( &device_data, need_lapl );
      lwd->eval_mmat_mgga_vxc_rks( &device_data, need_lapl );
    }
    else if( func.is_gga() ) lwd->eval_zmat_gga_vxc_rks( &device_data );
    else                     lwd->eval_zmat_lda_vxc_rks( &device_data );

    // Increment VXC (LT)
    lwd->inc_vxc( &device_data, func.is_mgga() );

  } // Loop over batches of batches 

  // Symmetrize VXC in device memory
  lwd->symmetrize_vxc( &device_data );

}



template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                       value_type* VXC, int64_t ldvxc, value_type* EXC, 
                       value_type *N_EL,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data ) {

  // Get integrate and keep data on device
  exc_vxc_local_work_( basis, P, ldp, task_begin, task_end, device_data );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  rt.device_backend()->master_queue_synchronize();

  // Receive XC terms from host
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
    device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXC, ldvxc );
  });

}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                                const value_type* Pz, int64_t ldpz,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {
  GauXC::util::unused(basis,Ps,ldps,Pz,ldpz,task_begin,task_end,device_data);
  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            value_type* VXCs, int64_t ldvxcs,
                            value_type* VXCz, int64_t ldvxcz, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {

  GauXC::util::unused(basis,Ps,ldps,Pz,ldpz,VXCs,ldvxcs,VXCz,ldvxcz,EXC,N_EL,task_begin,task_end,device_data);
  GAUXC_GENERIC_EXCEPTION("UKS NOT YET IMPLEMENTED FOR DEVICE");
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {
  GauXC::util::unused(basis,Ps,ldps,Pz,ldpz,Py,ldpy,Px,ldpx,task_begin,task_end,device_data);
  GAUXC_GENERIC_EXCEPTION("GKS NOT YET IMPLEMENTED FOR DEVICE");
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

  GauXC::util::unused(basis,Ps,ldps,Pz,ldpz,Py,ldpy,Px,ldpx,VXCs,ldvxcs,VXCz,ldvxcz,VXCy,ldvxcy,VXCx,ldvxcx,EXC,N_EL,task_begin,task_end,device_data);
  GAUXC_GENERIC_EXCEPTION("GKS NOT YET IMPLEMENTED FOR DEVICE");
}


}
}

