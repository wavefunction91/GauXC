#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include <stdexcept>
#include "device/xc_device_aos_data.hpp"
#include <fstream>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* EXC_GRAD ) { 
                 
  const auto& basis = this->load_balancer_->basis();

  // Check that P is sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");

  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [=](){ return lwd->create_device_data(); });

  const auto& mol = this->load_balancer_->molecule();
  const auto natoms = mol.size();
  if( this->reduction_driver_->takes_device_memory() ) {
    GAUXC_GENERIC_EXCEPTION("Device Reduction + EXC Grad NYI");
  } else {

    // Compute local contributions to EXC Gradient and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork", [&](){
      eval_exc_grad_local_work_( basis, P, ldp, EXC_GRAD, tasks.begin(),
        tasks.end(), *device_data_ptr );
    });

    this->timer_.time_op("XCIntegrator.ImbalanceWait",[&](){
      MPI_Barrier(this->load_balancer_->comm());
    });  

    this->timer_.time_op("XCIntegrator.Allreduce", [&](){
      this->reduction_driver_->allreduce_inplace( EXC_GRAD, 3*natoms, 
        ReductionOp::Sum );
    });

  }

}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_grad_local_work_( const basis_type& basis, 
    const value_type* P, int64_t ldp,
    host_task_iterator task_begin, host_task_iterator task_end,
    XCDeviceData& device_data ) {

  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  // Get basis map
  BasisSetMap basis_map(basis,mol);

  // Populate submat maps
  device_data.populate_submat_maps( basis.nbf(), task_begin, task_end, basis_map );

  // Sort tasks 
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( task_begin, task_end, task_comparator );




  // TODO: Refactor this into separate function
  auto& lb_state = this->load_balancer_->state();

  // Modify weights if need be
  if( not lb_state.modified_weights_are_stored ) {

  integrator_term_tracker enabled_terms;
  enabled_terms.weights = true;

  this->timer_.time_op("XCIntegrator.Weights", [&]() { 
    const auto natoms = mol.natoms();
    device_data.reset_allocations();
    device_data.allocate_static_data_weights( natoms );
    device_data.send_static_data_weights( mol, meta );

    // Processes batches in groups that saturadate available device memory
    auto task_it = task_begin;
    while( task_it != task_end ) {
      
      // Determine next task batch, send relevant data to device (weights only)
      auto task_batch_end = 
        device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

      // Apply partition weights 
      lwd->partition_weights( &device_data );
      
      // Copy back to host data
      device_data.copy_weights_to_tasks( task_it, task_batch_end );

      // Update iterator
      task_it = task_batch_end;

    } // End loop over batches

    // Signal that we don't need to do weights again
    lb_state.modified_weights_are_stored = true;
  });

  } // Mofify Weights

  // Do XC integration in task batches
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  const auto natoms  = mol.size();
  device_data.reset_allocations();
  device_data.allocate_static_data_exc_grad( nbf, nshells, natoms );
  device_data.send_static_data_density_basis( P, ldp, basis );

  // Zero integrands
  device_data.zero_exc_grad_integrands();

  // Processes batches in groups that saturadate available device memory
  integrator_term_tracker enabled_terms;
  enabled_terms.exc_grad = true;
  if( func.is_lda() )      enabled_terms.xc_approx = integrator_xc_approx::LDA; 
  else if( func.is_gga() ) enabled_terms.xc_approx = integrator_xc_approx::GGA; 
  else GAUXC_GENERIC_EXCEPTION("XC Approx NYI");

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXC Gradient only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Evaluate collocation
    if( func.is_gga() ) lwd->eval_collocation_hessian ( &device_data );
    else                lwd->eval_collocation_gradient( &device_data );

    // Evaluate X matrix
    const bool do_xmat_grad = func.is_gga();
    lwd->eval_xmat( &device_data, do_xmat_grad );

    // Evaluate U/V variables
    if( func.is_gga() ) lwd->eval_uvvar_gga( &device_data );
    else                lwd->eval_uvvar_lda( &device_data );

    // Evaluate XC functional (we need VXC for EXC Gradient)
    if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                lwd->eval_kern_exc_vxc_lda( func, &device_data );

    // Do scalar N_EL integration
    lwd->inc_nel( &device_data );

    // Increment EXC Gradient
    if( func.is_gga() ) lwd->inc_exc_grad_gga( &device_data );
    else                lwd->inc_exc_grad_lda( &device_data );

  } // Loop over batches of batches 

}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_grad_local_work_( const basis_type& basis, 
    const value_type* P, int64_t ldp, value_type* EXC_GRAD, 
    host_task_iterator task_begin, host_task_iterator task_end,
    XCDeviceData& device_data ) {

  // Compute XC gradient and keep data on the device
  eval_exc_grad_local_work_( basis, P, ldp, task_begin, task_end, device_data );

  // Receive XC gradient from host
  double N_EL;
  device_data.retrieve_exc_grad_integrands( EXC_GRAD, &N_EL );

  std::cout << N_EL << std::endl;
}


}
}
