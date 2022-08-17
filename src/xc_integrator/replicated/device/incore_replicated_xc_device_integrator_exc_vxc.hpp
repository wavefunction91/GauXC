#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/exceptions.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {


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
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [=](){ return lwd->create_device_data(); });


  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  if( this->reduction_driver_->takes_device_memory() ) {

    // If we can do reductions on the device (e.g. NCCL)
    // Don't communicate data back to the hot before reduction
    this->timer_.time_op("XCIntegrator.LocalWork", [&](){
      exc_vxc_local_work_( basis, P, ldp, tasks.begin(), tasks.end(), 
        *device_data_ptr);
    });

    this->timer_.time_op("XCIntegrator.ImbalanceWait",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  

    // Reduce results in device memory
    auto vxc_device = device_data_ptr->vxc_device_data();
    auto exc_device = device_data_ptr->exc_device_data();
    auto nel_device = device_data_ptr->nel_device_data();
    auto queue      = device_data_ptr->queue();
    this->timer_.time_op("XCIntegrator.Allreduce", [&](){
      this->reduction_driver_->allreduce_inplace( vxc_device, nbf*nbf, ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( exc_device, 1,       ReductionOp::Sum, queue );
      this->reduction_driver_->allreduce_inplace( nel_device, 1,       ReductionOp::Sum, queue );
    });

    // Retrieve data to host
    this->timer_.time_op("XCIntegrator.DeviceToHostCopy",[&](){
      device_data_ptr->retrieve_exc_vxc_integrands( EXC, &N_EL, VXC, ldvxc );
    });


  } else {

    // Compute local contributions to EXC/VXC and retrieve
    // data from device 
    this->timer_.time_op("XCIntegrator.LocalWork", [&](){
      exc_vxc_local_work_( basis, P, ldp, VXC, ldvxc, EXC, 
        &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
    });

    this->timer_.time_op("XCIntegrator.ImbalanceWait",[&](){
      MPI_Barrier(this->load_balancer_->runtime().comm());
    });  

    // Reduce Results in host mem
    this->timer_.time_op("XCIntegrator.Allreduce", [&](){
      this->reduction_driver_->allreduce_inplace( VXC, nbf*nbf, ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
      this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );
    });

  }

  std::cout << std::scientific << std::setprecision(12);
  std::cout << "NEL = " << N_EL << std::endl;
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


#if 0
  // Allocate static data on the stack
  const auto natoms  = mol.natoms();
  const auto nbf     = basis.nbf();
  const auto nshells = basis.nshells();
  //device_data.allocate_static_data( natoms, nbf, nshells );
  device_data.allocate_static_data_weights( natoms );
  device_data.allocate_static_data_exc_vxc( nbf, nshells );

  // Copy static data to device
  //device_data.send_static_data( P, ldp, basis, mol, meta );
  device_data.send_static_data_weights( mol, meta );
  device_data.send_static_data_density_basis( P, ldp, basis );

  // Zero integrands
  device_data.zero_exc_vxc_integrands();

  integrator_term_tracker enabled_terms;
  enabled_terms.weights = true;
  enabled_terms.exc_vxc = true;

  auto& lb_state = this->load_balancer_->state();

  // Processes batches in groups that saturadate available device memory
  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device
    auto task_batch_end = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Apply partition weights
    if( not lb_state.modified_weights_are_stored ) {
      lwd->partition_weights( &device_data );
      // Copy back to host data
      device_data.copy_weights_to_tasks( task_it, task_batch_end );
    }

    task_it = task_batch_end; // Update iterator


    // Evaluate collocation
    if( func.is_gga() ) lwd->eval_collocation_gradient( &device_data );
    else                lwd->eval_collocation( &device_data );

    // Evaluate X matrix
    lwd->eval_xmat( &device_data );

    // Evaluate U/V variables
    if( func.is_gga() ) lwd->eval_uvvar_gga( &device_data );
    else                lwd->eval_uvvar_lda( &device_data );

    // Evaluate XC functional
    if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                lwd->eval_kern_exc_vxc_lda( func, &device_data );

    // Do scalar EXC/N_EL integrations
    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );

    // Evaluate Z matrix
    if( func.is_gga() ) lwd->eval_zmat_gga_vxc( &device_data );
    else                lwd->eval_zmat_gga_vxc( &device_data );

    // Increment VXC (LT)
    lwd->inc_vxc( &device_data );

  } // Loop over batches of batches 

  lb_state.modified_weights_are_stored = true;

  // Symmetrize VXC in device memory
  lwd->symmetrize_vxc( &device_data );
#else


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
  if( func.is_lda() )      enabled_terms.xc_approx = integrator_xc_approx::LDA; 
  else if( func.is_gga() ) enabled_terms.xc_approx = integrator_xc_approx::GGA; 
  else GAUXC_GENERIC_EXCEPTION("XC Approx NYI");

  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device (EXC VXC only)
    task_it = 
      device_data.generate_buffers( enabled_terms, basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Evaluate collocation
    if( func.is_gga() ) lwd->eval_collocation_gradient( &device_data );
    else                lwd->eval_collocation( &device_data );

    // Evaluate X matrix
    lwd->eval_xmat( &device_data );

    // Evaluate U/V variables
    if( func.is_gga() ) lwd->eval_uvvar_gga( &device_data );
    else                lwd->eval_uvvar_lda( &device_data );

    // Evaluate XC functional
    if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                lwd->eval_kern_exc_vxc_lda( func, &device_data );

    // Do scalar EXC/N_EL integrations
    lwd->inc_exc( &device_data );
    lwd->inc_nel( &device_data );

    // Evaluate Z matrix
    if( func.is_gga() ) lwd->eval_zmat_gga_vxc( &device_data );
    else                lwd->eval_zmat_lda_vxc( &device_data );

    // Increment VXC (LT)
    lwd->inc_vxc( &device_data );

  } // Loop over batches of batches 

  // Symmetrize VXC in device memory
  lwd->symmetrize_vxc( &device_data );
#endif

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

  // Receive XC terms from host
  this->timer_.time_op("XCIntegrator.DeviceToHostCopy",[&](){
    device_data.retrieve_exc_vxc_integrands( EXC, N_EL, VXC, ldvxc );
  });

}

}
}

