/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "shellbatched_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"
#include "device/xc_device_aos_data.hpp"
#include <fstream>
#include <gauxc/exceptions.hpp>
#include <gauxc/util/unused.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_( int64_t m, int64_t n, const value_type* P, int64_t ldp, 
             const value_type* Pz, int64_t,
             const value_type*, int64_t,
             const value_type*, int64_t,
             value_type* EXC, const IntegratorSettingsXC& settings ) {


  if(Pz) GAUXC_GENERIC_EXCEPTION("UKS/GKS + EXC Only GPU NYI");

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");


  // Get Tasks
  auto& tasks = this->load_balancer_->get_tasks();

  // Allocate Device memory
  auto* lwd = dynamic_cast<LocalDeviceWorkDriver*>(this->local_work_driver_.get() );
  auto rt  = detail::as_device_runtime(this->load_balancer_->runtime());
  auto device_data_ptr = 
    this->timer_.time_op("XCIntegrator.DeviceAlloc",
      [&](){ return lwd->create_device_data(rt); });

  // Generate incore integrator instance, transfer ownership of LWD
  incore_integrator_type incore_integrator( this->func_, this->load_balancer_,
    this->release_local_work_driver(), this->reduction_driver_ );

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute local contributions to EXC/VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( basis, P, ldp, nullptr, 0, EXC, 
      &N_EL, tasks.begin(), tasks.end(), incore_integrator,
      *device_data_ptr );
  });

  // Release ownership of LWD back to this integrator instance
  this->local_work_driver_ = std::move( incore_integrator.release_local_work_driver() );


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );
  });


}



template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps, 
             const value_type* Pz, int64_t ldpz,
             value_type* EXC, const IntegratorSettingsXC& settings ) {

  eval_exc_(m, n, Ps, ldps, Pz, ldpz, nullptr, 0, nullptr, 0, EXC, settings);

}

template <typename ValueType>
void ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::
  eval_exc_( int64_t m, int64_t n, const value_type* P, int64_t ldp, 
             value_type* EXC, const IntegratorSettingsXC& settings ) {

  eval_exc_(m, n, P, ldp, nullptr, 0, nullptr, 0, nullptr, 0, EXC, settings);

}


}
}

