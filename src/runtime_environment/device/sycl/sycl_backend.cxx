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
#include "sycl_backend.hpp"

namespace GauXC {

SYCLBackend::SYCLBackend() {

  // Create SYCL queue. oneMKL is queue-driven, so the "BLAS handle" is a
  // binding to the same queue rather than a separate object to be tied to it
  device  = util::sycl_default_device();
  context = ::sycl::context( device );

  master_stream = std::make_shared< util::sycl_queue >( context, device );
  master_handle = std::make_shared< util::onemkl_handle >( master_stream );

#ifdef GAUXC_HAS_MAGMA
  #error "MAGMA is not supported with the SYCL backend"
#endif

}

SYCLBackend::~SYCLBackend() noexcept = default;

SYCLBackend::device_buffer_t SYCLBackend::allocate_device_buffer(int64_t sz) {
  void* ptr = nullptr;
  GAUXC_SYCL_ERROR( "SYCL Malloc Failed",
    ptr = ::sycl::malloc_device( sz, master_stream->queue ) );
  if( !ptr ) GAUXC_GENERIC_EXCEPTION("SYCL Malloc Returned NULL");
  return device_buffer_t{ptr,sz};
}

size_t SYCLBackend::get_available_mem() {
  // free_memory requires ZES_ENABLE_SYSMAN=1 on Level Zero; fall back to the
  // total global memory when the query is unavailable
  if( device.has( ::sycl::aspect::ext_intel_free_memory ) )
    return device.get_info< ::sycl::ext::intel::info::device::free_memory >();
  else
    return device.get_info< ::sycl::info::device::global_mem_size >();
}

void SYCLBackend::free_device_buffer( void* ptr ) {
  // cudaFree implicitly synchronizes the device, so the shared code is free
  // to release a buffer without first waiting on the work that reads it.
  // sycl::free carries no such guarantee -- freeing a USM allocation while a
  // kernel still references it is undefined behaviour -- so the outstanding
  // work is drained first to preserve the same contract.
  master_stream->synchronize();
  GAUXC_SYCL_ERROR( "Free Failed", ::sycl::free( ptr, master_stream->queue ) );
}

void SYCLBackend::master_queue_synchronize() {
  master_stream->synchronize();
}

device_queue SYCLBackend::queue() {
  return device_queue(master_stream);
}

void SYCLBackend::create_blas_queue_pool(int32_t ns) {
  blas_streams.resize(ns);
  blas_handles.resize(ns);
  for( auto i = 0; i < ns; ++i ) {
    blas_streams[i] = std::make_shared<util::sycl_queue>( context, device );
    blas_handles[i] = std::make_shared<util::onemkl_handle>( blas_streams[i] );
  }
}

void SYCLBackend::sync_master_with_blas_pool() {
  const auto n_streams = blas_streams.size();
  std::vector<util::sycl_event> blas_events( n_streams );
  for( size_t iS = 0; iS < n_streams; ++iS )
    blas_events[iS].record( *blas_streams[iS] );

  for( auto& event : blas_events ) master_stream->wait(event);
}

void SYCLBackend::sync_blas_pool_with_master() {
  util::sycl_event master_event;
  master_event.record( *master_stream );
  for( auto& stream : blas_streams ) stream->wait( master_event );
}

size_t SYCLBackend::blas_pool_size(){ return blas_streams.size(); }

device_queue SYCLBackend::blas_pool_queue(int32_t i) {
  return device_queue( blas_streams.at(i) );
}

device_blas_handle SYCLBackend::blas_pool_handle(int32_t i) {
  return device_blas_handle( blas_handles.at(i) );
}
device_blas_handle SYCLBackend::master_blas_handle() {
  return device_blas_handle( master_handle );
}

void SYCLBackend::copy_async_( size_t sz, const void* src, void* dest,
  std::string msg ) {
  // Despite the name, this must complete before returning.
  //
  // cudaMemcpyAsync is only truly asynchronous when the host side is pinned
  // memory; for a device-to-pageable-host transfer CUDA makes it implicitly
  // synchronizing (the driver stages through an internal buffer and blocks).
  // The shared code relies on that: retrieve_exc_vxc_integrands issues these
  // copies into ordinary host allocations and its caller reads them straight
  // away, with no synchronization in between.
  //
  // SYCL's queue::memcpy offers no such guarantee -- it returns an event and
  // nothing more -- so without waiting here the host races the in-flight copy
  // and reads garbage (observed as NaN in VXC, then heap corruption).
  GAUXC_SYCL_ERROR( "SYCL Memcpy Async Failed ["+msg+"]",
    master_stream->queue.memcpy( dest, src, sz ).wait_and_throw() );
}

void SYCLBackend::set_zero_(size_t sz, void* data, std::string msg ) {
  GAUXC_SYCL_ERROR( "SYCL Memset Failed ["+msg+"]",
    master_stream->queue.memset( data, 0, sz ).wait_and_throw() );
}

void SYCLBackend::set_zero_async_master_queue_(size_t sz, void* data, std::string msg ) {
  GAUXC_SYCL_ERROR( "SYCL Memset Failed ["+msg+"]",
    master_stream->queue.memset( data, 0, sz ) );
}

void SYCLBackend::copy_async_2d_( size_t M, size_t N, const void* A, size_t LDA,
  void* B, size_t LDB, std::string msg ) {

  // A strided copy of N rows, M bytes each.
  //
  // queue::ext_oneapi_memcpy2d is deliberately avoided: its helper kernel is
  // not emitted into an ahead-of-time device image, so an AoT build faults
  // inside the runtime's kernel-info lookup. The same entry point is also
  // implicated in a known Level Zero defect on this platform
  // (CMPLRLLVM-76543 / GSD-11506). One enqueue per row keeps the copy on the
  // ordinary USM memcpy path, which both AoT and JIT handle.
  //
  // The copy is completed before returning, for the same reason as
  // copy_async_: callers inherit CUDA's implicitly-synchronizing semantics
  // for a device-to-pageable-host transfer.
  if( !M or !N ) return;

  auto* dst = static_cast<unsigned char*>(B);
  const auto* src = static_cast<const unsigned char*>(A);

  GAUXC_SYCL_ERROR( "SYCL 2D Memcpy Async Failed ["+msg+"]",
    [&]() {
      // Contiguous on both sides: a single memcpy covers the whole region.
      if( LDA == M and LDB == M ) {
        master_stream->queue.memcpy( dst, src, M*N );
      } else {
        for( size_t row = 0; row < N; ++row )
          master_stream->queue.memcpy( dst + row*LDB, src + row*LDA, M );
      }
      master_stream->queue.wait_and_throw();
    }()
  );
}


void SYCLBackend::check_error_(std::string msg) {
  // SYCL reports kernel failures asynchronously through the queue's exception
  // handler. Draining it here is the closest analogue to cudaGetLastError
  GAUXC_SYCL_ERROR( "SYCL Failed ["+msg+"]",
    master_stream->queue.throw_asynchronous() );
}

std::unique_ptr<DeviceBackend> make_device_backend() {
  return std::make_unique<SYCLBackend>();
}
}
