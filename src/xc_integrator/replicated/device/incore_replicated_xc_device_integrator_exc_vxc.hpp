/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
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
      //exc_vxc_local_work_( basis, P, ldp, tasks.begin(), tasks.end(), *device_data_ptr);
      exc_vxc_local_work_( basis, P, ldp, nullptr, 0, nullptr, 0, nullptr, 0, tasks.begin(), tasks.end(), *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce results in device memory
    auto vxc_device = device_data_ptr->vxc_s_device_data();
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
      exc_vxc_local_work_( basis, P, ldp, nullptr, 0, nullptr, 0, nullptr, 0,
                                  VXC, ldvxc, nullptr, 0, nullptr, 0, nullptr, 0,
                                  EXC, &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
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
  //GauXC::util::unused(m,n,Ps,ldps,Pz,ldpz,VXCs,ldvxcs,VXCz,ldvxcz,EXC,settings);

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
  if( ldpz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPz");
  if( ldvxcz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCz");


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
      exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, tasks.begin(), tasks.end(), 
        *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce results in device memory
    auto vxc_s_device = device_data_ptr->vxc_s_device_data();
    auto vxc_z_device = device_data_ptr->vxc_z_device_data();
    auto exc_device = device_data_ptr->exc_device_data();
    auto nel_device = device_data_ptr->nel_device_data();
    auto queue      = device_data_ptr->queue();
    this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
      this->reduction_driver_->allreduce_inplace( vxc_s_device, nbf*nbf, ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( vxc_z_device, nbf*nbf, ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
    });

    // Retrieve data to host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
      device_data_ptr->retrieve_exc_vxc_integrands( EXC, &N_EL, VXCs, ldvxcs, VXCz, ldvxcz );
    });


  } else {

    // Compute local contributions to EXC/VXC and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork_EXC_VXC", [&](){
      exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0,
                                VXCs, ldvxcs, VXCz, ldvxcz, nullptr, 0, nullptr, 0, EXC, 
                              &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
    });

    GAUXC_MPI_CODE(
    this->timer_.time_op("XCIntegrator.ImbalanceWait_EXC_VXC",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  
    )

    // Reduce Results in host mem
    this->timer_.time_op("XCIntegrator.Allreduce_EXC_VXC", [&](){
      this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
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
  exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data ) {
  const bool is_gks = (Pz != nullptr) and (Py != nullptr) and (Px != nullptr);
  const bool is_uks = (Pz != nullptr) and (Py == nullptr) and (Px == nullptr);
  const bool is_rks = not is_uks and not is_gks;
  if (not is_rks and not is_uks and not is_gks) {
    GAUXC_GENERIC_EXCEPTION("MUST BE EITHER RKS, UKS, or GKS!");
  }
  
  if (is_gks) GAUXC_GENERIC_EXCEPTION( "GKS DEVICE NYI!");

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
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified"); 
  }
  

  integrator_term_tracker enabled_terms;
  enabled_terms.exc_vxc = true;
  if (is_rks) enabled_terms.ks_scheme = RKS;
  if (is_uks) enabled_terms.ks_scheme = UKS;
  if (is_gks) enabled_terms.ks_scheme = GKS;
  
  // Do XC integration in task batches
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  device_data.reset_allocations();
  device_data.allocate_static_data_exc_vxc( nbf, nshells, enabled_terms );
  
  if (is_rks) device_data.send_static_data_density_basis( Ps, ldps, basis );
  else if (is_uks) device_data.send_static_data_density_basis( Ps, ldps, Pz, ldpz, basis );
  //if (is_gks) device_data.send_static_data_density_basis( Ps, ldps, Pz, ldpz, Px, ldpx, Py, ldpy, basis );

    // for debugging
    auto* data = dynamic_cast<XCDeviceStackData*>(&device_data);
    auto base_stack = data->base_stack;
    auto static_stack = data->static_stack;

  // Processes batches in groups that saturate available device memory
  if( func.is_lda() )      enabled_terms.xc_approx = integrator_xc_approx::LDA; 
  else if( func.is_gga() ) enabled_terms.xc_approx = integrator_xc_approx::GGA; 
  else GAUXC_GENERIC_EXCEPTION("XC Approx NYI");

  // Zero integrands
  device_data.zero_exc_vxc_integrands(enabled_terms);
  

  if( func.is_gga() and is_uks ) GAUXC_GENERIC_EXCEPTION("UKS GGA NYI");

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXC VXC only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/
    
    // Evaluate collocation
    if( func.is_gga() ) lwd->eval_collocation_gradient( &device_data );
    else                lwd->eval_collocation( &device_data );

    
    double xmat_fac = 1.0;
    if (is_rks) {
      xmat_fac = 2.0;
      // Evaluate X matrix
      lwd->eval_xmat( xmat_fac, &device_data, false, DEN_S );
      
    }

    else if (is_uks) {
      xmat_fac = 0.5;
      // Evaluate X matrix
      lwd->eval_xmat( xmat_fac, &device_data, false, DEN_S );
      // Contract X matrix with bf -> den_eval
      lwd->eval_den( &device_data, DEN_S );
      // Repeat for Z density
      lwd->eval_xmat( xmat_fac, &device_data, false, DEN_Z );
      lwd->eval_den( &device_data, DEN_Z );

      // Evaluate U/V variables
      if( func.is_gga() ) GAUXC_GENERIC_EXCEPTION("Device UKS+GGA NYI!");
    }

    // Evaluate U/V variables
    if( func.is_gga() ) lwd->eval_uvvar_gga( &device_data, enabled_terms );
    else                lwd->eval_uvvar_lda( &device_data, enabled_terms );

    // Evaluate XC functional
    if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                lwd->eval_kern_exc_vxc_lda( func, &device_data );
    

    // Do scalar EXC/N_EL integrations
    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );


    if (is_rks) {
      // Evaluate Z matrix
      if( func.is_gga() ) lwd->eval_zmat_gga_vxc_rks( &device_data );
      else                lwd->eval_zmat_lda_vxc_rks( &device_data );
      // Increment VXC
      lwd->inc_vxc( &device_data, DEN_S );
    }
    if (is_uks) {
      // Evaluate Scalar Z matrix
      //if( func.is_gga() ) lwd->eval_zmat_gga_vxc_uks( &device_data, DEN_S );
      if( func.is_gga() ) GAUXC_GENERIC_EXCEPTION("UKS GGA eval_zmat NYI");
      else                lwd->eval_zmat_lda_vxc_uks( &device_data, DEN_S );
      // Increment Scalar VXC
      lwd->inc_vxc( &device_data, DEN_S );
      // Repeat for Z VXC

      //if( func.is_gga() ) lwd->eval_zmat_gga_vxc_uks( &device_data, DEN_Z );
      if( func.is_gga() ) GAUXC_GENERIC_EXCEPTION("UKS GGA eval_zmat NYI");
      else                lwd->eval_zmat_lda_vxc_uks( &device_data, DEN_Z );
      lwd->inc_vxc( &device_data, DEN_Z );
      
    }

  } // Loop over batches of batches 

  // Symmetrize VXC in device memory

  lwd->symmetrize_vxc( &device_data, DEN_S );
  if (is_uks) {
    lwd->symmetrize_vxc( &device_data, DEN_Z );
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
  
  //GauXC::util::unused(basis,Ps,ldps,Pz,ldpz,Py,ldpy,Px,ldpx,VXCs,ldvxcs,VXCz,ldvxcz,VXCy,ldvxcy,VXCx,ldvxcx,EXC,N_EL,task_begin,task_end,device_data);
  const bool is_gks = (Pz != nullptr) and (Py != nullptr) and (Px != nullptr);
  const bool is_uks = (Pz != nullptr) and (Py == nullptr) and (Px == nullptr);
  const bool is_rks = not is_uks and not is_gks;
  if (not is_rks and not is_uks and not is_gks) {
    GAUXC_GENERIC_EXCEPTION("MUST BE EITHER RKS, UKS, or GKS!");
  }
  

  // Get integrate and keep data on device
  exc_vxc_local_work_( basis, Ps, ldps, Pz, ldpz, Py, ldpy, Px, ldpx, task_begin, task_end, device_data );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  rt.device_backend()->master_queue_synchronize();

  // Receive XC terms from host
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy_EXC_VXC",[&](){
    if (is_rks)  device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXCs, ldvxcs ); 
    if (is_uks)  device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXCs, ldvxcs, VXCz, ldvxcz ); 
    //if (is_gks)  device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXCs, ldvxcs, VXCz, ldvxcz, VXCy, ldvxcy, VXCx, ldvxcx ); 
  });

}


}
}

