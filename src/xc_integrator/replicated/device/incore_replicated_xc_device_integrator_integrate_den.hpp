#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/exceptions.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  integrate_den_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* N_EL ) {


  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
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
  auto& rt  = dynamic_cast<const DeviceRuntimeEnvironment&>(this->load_balancer_->runtime());
  auto device_data_ptr = lwd->create_device_data(rt);


  if( this->reduction_driver_->takes_device_memory() ) {
    GAUXC_GENERIC_EXCEPTION("Device Reduction + Integrate Den NYI");
  } else {

    // Compute local contributions to N_EL and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork_IntegrateDen", [&](){
      integrate_den_local_work_( basis, P, ldp, N_EL,
        tasks.begin(), tasks.end(), *device_data_ptr);
    });

    this->timer_.time_op("XCIntegrator.ImbalanceWait",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  

    // Reduce Results in host mem
    this->timer_.time_op("XCIntegrator.Allreduce", [&](){
      this->reduction_driver_->allreduce_inplace( N_EL, 1, ReductionOp::Sum );
    });

  }
}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  integrate_den_local_work_( const basis_type& basis, const value_type* P, 
                       int64_t ldp, value_type* N_EL,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data ) {
#if 1
  GAUXC_GENERIC_EXCEPTION("NYI");
#else


  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

  // Setup Aliases
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

  }

  // Do XC integration in task batches
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  device_data.reset_allocations();
  device_data.allocate_static_data_integrate_den( nbf, nshells );
  device_data.send_static_data_density_basis( P, ldp, basis );

  // Zero integrands
  device_data.zero_integrate_den_integrands();

  // Processes batches in groups that saturadate available device memory
  integrator_term_tracker enabled_terms;
  enabled_terms.integrate_den = true;

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (Density only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Evaluate collocation
    lwd->eval_collocation( &device_data );

    // Evaluate X matrix
    lwd->eval_xmat( &device_data );

    // Do scalar N_EL integration
    lwd->inc_nel( &device_data );

  } // Loop over batches of batches 

  // Receive N_EL from device
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy",[&](){
    device_data.retrieve_integrate_den_integrands( N_EL );
  });
#endif
}

}
}

