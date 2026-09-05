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
#include "oneccl_reduction_driver.hpp"
#include <cstring>
#include <memory>
#include <map>
#include <mutex>
#include <iostream>

#include "device/sycl/sycl_backend.hpp"
#include "device/device_queue.hpp"

namespace GauXC {

ccl::datatype get_ccl_datatype( std::type_index idx ) {

  static std::map<std::type_index, ccl::datatype> map {
    {std::type_index(typeid(double)), ccl::datatype::float64},
    {std::type_index(typeid(float)),  ccl::datatype::float32}
  };

  return map.at(idx);

}

ccl::reduction get_ccl_op( ReductionOp op ) {

  static std::map< ReductionOp, ccl::reduction > map {
    { ReductionOp::Sum, ccl::reduction::sum }
  };

  return map.at(op);

}

// Returns the SYCLBackend bound to this runtime, or throws. oneCCL's
// create_communicator (unlike ncclCommInitRank) must be bound to a concrete
// sycl::device / sycl::context at construction time, so the reduction driver
// needs backend access up front rather than only at collective-call time.
SYCLBackend* get_sycl_backend_from_runtime( const RuntimeEnvironment& rt ) {
  auto* backend =
    dynamic_cast<SYCLBackend*>(detail::as_device_runtime(rt).device_backend());
  if( !backend ) GAUXC_GENERIC_EXCEPTION("OneCCL Reduction Driver Requires A SYCL Backend");
  return backend;
}

// SYCL equivalent of get_cuda_stream_from_optional_args: pulls a
// util::sycl_queue out of the type-erased device_queue when the caller
// supplied one. The CUDA driver's "stream == 0 -> cudaDeviceSynchronize()"
// fallback has no SYCL analogue (there is no default/null queue), so when no
// queue is supplied we fall back to the SYCLBackend's own master queue -
// the same queue whose device/context were used to build the communicator -
// rather than synchronizing the whole device.
util::sycl_queue* get_sycl_queue_from_optional_args( std::any& args,
  SYCLBackend* backend ) {

  if( args.has_value() ) {
    if( auto ptr = std::any_cast<device_queue>( &args ) )
    if( auto passed_stream = ptr->queue_as_ptr<util::sycl_queue>() ) {
      return passed_stream;
    }
  }

  device_queue backend_queue = backend->queue();
  auto* q = backend_queue.queue_as_ptr<util::sycl_queue>();
  if( !q ) GAUXC_GENERIC_EXCEPTION("SYCL Backend Did Not Return A Valid Queue");
  return q;
}


OneCCLReductionDriver::OneCCLReductionDriver(const RuntimeEnvironment& rt) :
  DeviceReductionDriver(rt) {

  auto* backend = get_sycl_backend_from_runtime(rt);
  oneccl_comm_ = std::make_shared<util::oneccl_comm>( rt.comm(),
    backend->device, backend->context );

}


OneCCLReductionDriver::~OneCCLReductionDriver() noexcept = default;
OneCCLReductionDriver::OneCCLReductionDriver(const OneCCLReductionDriver&) = default;


void OneCCLReductionDriver::allreduce_typeerased( const void* src, void* dest,
  size_t size, ReductionOp op, std::type_index idx, std::any optional_args )  {

  auto* backend = get_sycl_backend_from_runtime( runtime_ );
  auto* queue   = get_sycl_queue_from_optional_args( optional_args, backend );

  auto synchronize = [&]() { queue->synchronize(); };

  synchronize();
  try {
    auto stream = ccl::create_stream( queue->queue );
    ccl::allreduce( src, dest, size, get_ccl_datatype(idx),
      get_ccl_op(op), *oneccl_comm_, stream ).wait();
  } catch( const ccl::exception& e ) {
    GAUXC_GENERIC_EXCEPTION( std::string("OneCCL FAILED: ") + e.what() );
  }
  synchronize();

}
void OneCCLReductionDriver::allreduce_inplace_typeerased( void* data, size_t size,
  ReductionOp op, std::type_index idx, std::any optional_args) {

  auto* backend = get_sycl_backend_from_runtime( runtime_ );
  auto* queue   = get_sycl_queue_from_optional_args( optional_args, backend );

  auto synchronize = [&]() { queue->synchronize(); };

  synchronize();
  try {
    auto stream = ccl::create_stream( queue->queue );
    ccl::allreduce( data, data, size, get_ccl_datatype(idx),
      get_ccl_op(op), *oneccl_comm_, stream ).wait();
  } catch( const ccl::exception& e ) {
    GAUXC_GENERIC_EXCEPTION( std::string("OneCCL FAILED: ") + e.what() );
  }
  synchronize();

}

std::unique_ptr<detail::ReductionDriverImpl> OneCCLReductionDriver::clone() {
  return std::make_unique<OneCCLReductionDriver>(*this);
}


}
