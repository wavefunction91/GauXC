#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include <stdexcept>

namespace GauXC  {
namespace detail {

template <typename ValueType>
IncoreReplicatedXCDeviceIntegrator<ValueType>::~IncoreReplicatedXCDeviceIntegrator() noexcept = default;


template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {


  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  std::string fun_name = __PRETTY_FUNCTION__;
  if( m != n ) 
    throw std::logic_error(fun_name + " P/VXC Must Be Square");
  if( m != nbf ) 
    throw std::logic_error(fun_name + " P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    throw std::logic_error(fun_name + " Invalid LDP");
  if( ldvxc < nbf )
    throw std::logic_error(fun_name + " Invalid LDVXC");


  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

  // TODO: This is a bad name for this
  size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  // Allocate Device memory
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [=](){ return this->create_device_data(); });


  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute local contributions to EXC/VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( basis, P, ldp, VXC, ldvxc, EXC, 
      &N_EL, tasks.begin(), tasks.end(), *device_data_ptr);
  });

#ifdef GAUXC_ENABLE_MPI

  int world_size;
  auto comm = this->load_balancer_->comm();
  MPI_Comm_size( comm, &world_size );

  if( world_size > 1 ) {

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    // Test of communicator is an inter-communicator
    // XXX: Can't think of a case when this would be true, but who knows...
    int inter_flag;
    MPI_Comm_test_inter( comm, &inter_flag );

    // Is Intra-communicator, Allreduce can be done inplace
    if( not inter_flag ) {

      MPI_Allreduce( MPI_IN_PLACE, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, comm );
      MPI_Allreduce( MPI_IN_PLACE, EXC,  1, MPI_DOUBLE, MPI_SUM, comm );
      MPI_Allreduce( MPI_IN_PLACE, &N_EL, 1, MPI_DOUBLE, MPI_SUM, comm );

    // Isn't Intra-communicator (weird), Allreduce can't be done inplace
    } else {

      std::allocator<value_type> alloc;
      auto VXC_cpy = alloc.allocate( nbf*nbf );
      value_type EXC_cpy = *EXC, N_EL_cpy = N_EL;

      MPI_Allreduce( VXC_cpy, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, comm );
      MPI_Allreduce( &EXC_cpy,  EXC,  1, MPI_DOUBLE, MPI_SUM, comm );
      MPI_Allreduce( &N_EL_cpy, &N_EL, 1, MPI_DOUBLE, MPI_SUM, comm );
      

    }
  });

  }

#endif

}

template <typename ValueType>
void IncoreReplicatedXCDeviceIntegrator<ValueType>::
  exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                       value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                       host_task_iterator task_begin, host_task_iterator task_end,
                       XCDeviceData& device_data ) {


  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& mol   = this->load_balancer_->molecule();
  const auto& meta  = this->load_balancer_->molmeta();

  // Get basis map
  BasisSetMap basis_map(basis);

  // Sort tasks 
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.nbe) > (b.points.size() * b.nbe);
  };
  std::sort( task_begin, task_end, task_comparator );

  // TODO: This is a bad name for this
  const size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  // Allocate static data on the stack
  const auto natoms  = mol.natoms();
  const auto nbf     = basis.nbf();
  const auto nshells = basis.size();
  device_data.allocate_static_data( natoms, n_deriv, nbf, nshells );

  // TODO: Copy static data to device


  // Processes batches in groups that saturadate available device memory
  auto task_it = task_begin;
  while( task_it != task_end ) {

    // Determine next task batch, send relevant data to device
    task_it = device_data.generate_buffers( basis_map, task_it, task_end );

    /*** Process the batches ***/

    // Apply partition weights
    lwd->partition_weights( &device_data );

    // Evaluate collocation
    if( n_deriv == 1 ) lwd->eval_collocation_gradient( &device_data );
    else               lwd->eval_collocation( &device_data );

    // Evaluate X matrix
    lwd->eval_xmat( &device_data );

    // Evaluate U/V variables
    if( func.is_gga() ) lwd->eval_uvvar_gga( &device_data );
    else                lwd->eval_uvvar_lda( &device_data );

    // Evaluate XC functional
    if( func.is_gga() ) lwd->eval_kern_exc_vxc_gga( func, &device_data );
    else                lwd->eval_kern_exc_vxc_lda( func, &device_data );

    // TODO: Scalar integrations

    // Evaluate Z matrix
    if( func.is_gga() ) lwd->eval_zmat_gga_vxc( &device_data );
    else                lwd->eval_zmat_gga_vxc( &device_data );

    // Increment VXC (LT)
    lwd->inc_vxc( &device_data );

  } // Loop over batches of batches 


  // TODO: Receive XC terms from host
  // TODO: Symmetrize VXC
}


}
}
