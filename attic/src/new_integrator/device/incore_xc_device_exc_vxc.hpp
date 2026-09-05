#include <memory>
#include <gauxc/new_xc_integrator/replicated/incore_xc_device_integrator.hpp>

#include "device/local_work_replicated_incore_exc_vxc.hpp"
#include <gauxc/util/unused.hpp>
#include <gauxc/gauxc_config.hpp>
#include "exceptions/magma_exception.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
void IncoreXCDeviceIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {

#ifdef GAUXC_ENABLE_MAGMA
  // Initialize MAGMA
  {
    auto ierr = magma_init();
    GAUXC_MAGMA_ERROR( "MAGMA Init Failed", ierr );
  }
#endif

  util::unused(m,n,ldp,ldvxc);

  size_t nbf = this->basis_->nbf();

  //// TODO: Check that P is sane


  // Generate Tasks
  auto& tasks = this->load_balancer_->get_tasks();
  size_t n_deriv = this->func_->is_gga() ? 1 : 0;

  // Allocate Memory
  auto device_data = this->timer_.time_op("XCIntegrator.DeviceAlloc",[&]() {
    auto ptr = GauXC::integrator::device::make_device_data<ValueType>();
    ptr->allocate_static_data( this->load_balancer_->molecule().natoms(),
                               n_deriv, this->basis_->nbf(), 
                               this->basis_->size() );
    return ptr;
  });

  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    GauXC::integrator::device::local_work_replicated_incore_exc_vxc< value_type >(
      n_deriv, XCWeightAlg::SSF, state_, *this->func_, 
      *this->basis_, this->load_balancer_->molecule(), 
      this->load_balancer_->molmeta(), *device_data, tasks, P, 
      VXC, EXC, &N_EL 
    );
  });

  // Update State of Integrator
  state_.load_balancer_populated     = true;
  //state_.modified_weights_are_stored = true;

            
#ifdef GAUXC_ENABLE_MPI

  int world_size;
  MPI_Comm_size( this->comm_, &world_size );

  if( world_size > 1 ) {

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    // Test of communicator is an inter-communicator
    // XXX: Can't think of a case when this would be true, but who knows...
    int inter_flag;
    MPI_Comm_test_inter( this->comm_, &inter_flag );

    // Is Intra-communicator, Allreduce can be done inplace
    if( not inter_flag ) {

      MPI_Allreduce( MPI_IN_PLACE, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, this->comm_ );
      MPI_Allreduce( MPI_IN_PLACE, EXC,  1, MPI_DOUBLE, MPI_SUM, this->comm_ );
      MPI_Allreduce( MPI_IN_PLACE, &N_EL, 1, MPI_DOUBLE, MPI_SUM, this->comm_ );

    // Isn't Intra-communicator (weird), Allreduce can't be done inplace
    } else {

      std::allocator<value_type> alloc;
      auto VXC_cpy = alloc.allocate( nbf*nbf );
      value_type EXC_cpy = *EXC, N_EL_cpy = N_EL;

      MPI_Allreduce( VXC_cpy, VXC, nbf*nbf, MPI_DOUBLE,
                     MPI_SUM, this->comm_ );
      MPI_Allreduce( &EXC_cpy,  EXC,  1, MPI_DOUBLE, MPI_SUM, this->comm_ );
      MPI_Allreduce( &N_EL_cpy, &N_EL, 1, MPI_DOUBLE, MPI_SUM, this->comm_ );
      

    }
  });

  }

#endif

#ifdef GAUXC_ENABLE_MAGMA
  // Finalize MAGMA
  {
    auto ierr = magma_finalize();
    GAUXC_MAGMA_ERROR( "MAGMA Finalize Failed", ierr );
  }
#endif

}

}
}

