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
#pragma once
#include <gauxc/gauxc_config.hpp>
#include "exceptions/sycl_exception.hpp"

#ifdef GAUXC_HAS_SYCL
#include <sycl/sycl.hpp>
#include <utility>

namespace GauXC {
namespace util  {

struct sycl_queue;
struct sycl_event;

/// Async exception handler which promotes SYCL errors to GauXC exceptions
inline void sycl_async_error_handler( ::sycl::exception_list exceptions ) {
  for( auto& e : exceptions ) {
    try { std::rethrow_exception(e); }
    catch( const ::sycl::exception& err ) {
      throw sycl_exception( __FILE__, __LINE__, "SYCL Async Error", err );
    }
  }
}

/// Device selection honoring ZE_AFFINITY_MASK / ONEAPI_DEVICE_SELECTOR. Falls
/// back to the default selector when no GPU is visible so that CPU / FPGA
/// backends remain usable.
inline ::sycl::device sycl_default_device() {
  try {
    return ::sycl::device( ::sycl::gpu_selector_v );
  } catch( const ::sycl::exception& ) {
    return ::sycl::device( ::sycl::default_selector_v );
  }
}

/// RAII wrapper for an in-order SYCL queue. In-order semantics give the same
/// serialized execution the CUDA/HIP backends get from a single stream.
struct sycl_queue {

  ::sycl::queue queue;

  inline sycl_queue() :
    sycl_queue( ::sycl::context(sycl_default_device()), sycl_default_device() ) { }

  inline sycl_queue( const ::sycl::context& ctx, const ::sycl::device& dev ) :
    queue( ctx, dev, sycl_async_error_handler,
           ::sycl::property_list{ ::sycl::property::queue::in_order{} } ) { }

  inline ~sycl_queue() noexcept = default;

  sycl_queue( const sycl_queue& ) = delete;
  sycl_queue( sycl_queue&& ) noexcept = default;

  inline operator ::sycl::queue&() { return queue; }
  inline operator const ::sycl::queue&() const { return queue; }

  inline ::sycl::queue* operator->() { return &queue; }

  /// Enqueue a barrier which blocks this queue until "event" has completed
  inline void wait( ::sycl::event event ) {
    GAUXC_SYCL_ERROR( "QUEUE WAIT FAILED",
      queue.ext_oneapi_submit_barrier({event}) );
  }

  inline void synchronize() {
    GAUXC_SYCL_ERROR( "SYCL Queue Synchronize Failed", queue.wait_and_throw() );
  }

};


/// RAII wrapper for a SYCL event. Recorded via a queue barrier which mirrors
/// the cudaEventRecord / hipEventRecord semantics.
struct sycl_event {

  ::sycl::event event;

  inline sycl_event() = default;
  inline ~sycl_event() noexcept = default;

  sycl_event( const sycl_event& ) = default;
  sycl_event( sycl_event&& ) noexcept = default;

  inline operator ::sycl::event() const { return event; }

  inline void record( ::sycl::queue& queue ) {
    GAUXC_SYCL_ERROR( "Event Record Failed",
      event = queue.ext_oneapi_submit_barrier() );
  }

  inline void record( sycl_queue& queue ) { record( queue.queue ); }

};





template <typename T>
inline T* sycl_malloc( size_t n, ::sycl::queue& queue ) {

  T* ptr = nullptr;
  GAUXC_SYCL_ERROR( "SYCL Malloc Failed",
    ptr = ::sycl::malloc_device<T>( n, queue ) );
  if( !ptr ) GAUXC_GENERIC_EXCEPTION("SYCL Malloc Returned NULL");

  return ptr;
}

template <typename T>
inline T* sycl_malloc_host( size_t n, ::sycl::queue& queue ) {

  T* ptr = nullptr;
  GAUXC_SYCL_ERROR( "SYCL Malloc Host Failed",
    ptr = ::sycl::malloc_host<T>( n, queue ) );
  if( !ptr ) GAUXC_GENERIC_EXCEPTION("SYCL Malloc Host Returned NULL");

  return ptr;
}

/// Shared (host + device accessible) USM allocation. Needed for parameter
/// arrays that a host-side dispatch routine reads (e.g. oneMKL's variable-group
/// gemm_batch, which enumerates group boundaries before enqueueing), as
/// opposed to device-only scratch that only device kernels touch.
template <typename T>
inline T* sycl_malloc_shared( size_t n, ::sycl::queue& queue ) {

  T* ptr = nullptr;
  GAUXC_SYCL_ERROR( "SYCL Malloc Shared Failed",
    ptr = ::sycl::malloc_shared<T>( n, queue ) );
  if( !ptr ) GAUXC_GENERIC_EXCEPTION("SYCL Malloc Shared Returned NULL");

  return ptr;
}






template <typename T>
inline void sycl_free( ::sycl::queue& queue, T*& ptr ) {
  // cuda_free's cudaFree implicitly synchronizes the device; sycl::free does
  // not, and releasing a USM allocation that enqueued work still references
  // is undefined behaviour. Drain the queue first so callers written against
  // the CUDA helper stay correct.
  GAUXC_SYCL_ERROR( "SYCL Free Sync Failed", queue.wait_and_throw() );
  GAUXC_SYCL_ERROR( "SYCL Free Failed", ::sycl::free( (void*)ptr, queue ) );
  ptr = nullptr;
}

template <typename T, typename... Args>
inline void sycl_free( ::sycl::queue& queue, T*& ptr, Args&&... args ) {
  sycl_free(queue, ptr);
  sycl_free(queue, std::forward<Args>(args)...);
}




template <typename T>
inline void sycl_copy( size_t len, T* dest, const T* src, ::sycl::queue& queue,
                       std::string m = "") {
  GAUXC_SYCL_ERROR( "SYCL Memcpy Failed ["+m+"]",
    queue.memcpy( dest, src, len * sizeof(T) ).wait_and_throw() );
}

template <typename T>
inline void sycl_copy_async( size_t len, T* dest, const T* src,
                             ::sycl::queue& queue, std::string m = "" ) {
  GAUXC_SYCL_ERROR( "SYCL Memcpy Async Failed ["+m+"]",
    queue.memcpy( dest, src, len * sizeof(T) ) );
}


template <typename T>
inline void sycl_copy_2d( T* dest, size_t dest_pitch, const T* src,
                          size_t src_pitch, size_t width, size_t height,
                          ::sycl::queue& queue, std::string m = "" ) {
  GAUXC_SYCL_ERROR( "SYCL 2D Memcpy Failed ["+m+"]",
    queue.ext_oneapi_memcpy2d( dest, dest_pitch, src, src_pitch, width,
      height ).wait_and_throw() );
}

template <typename T>
inline void sycl_copy_2d_async( T* dest, size_t dest_pitch, const T* src,
                                size_t src_pitch, size_t width, size_t height,
                                ::sycl::queue& queue, std::string m = "" ) {
  GAUXC_SYCL_ERROR( "SYCL 2D Memcpy Async Failed ["+m+"]",
    queue.ext_oneapi_memcpy2d( dest, dest_pitch, src, src_pitch, width, height ) );
}


template <typename T>
inline void sycl_set_zero( size_t len, T* ptr, ::sycl::queue& queue,
                           std::string m = "" ) {
  GAUXC_SYCL_ERROR( "SYCL Memset Failed ["+m+"]",
    queue.memset( ptr, 0, len * sizeof(T) ).wait_and_throw() );
}

template <typename T>
inline void sycl_set_zero_async( size_t len, T* ptr, ::sycl::queue& queue,
                                 std::string m = "" ) {
  GAUXC_SYCL_ERROR( "SYCL Memset Async Failed ["+m+"]",
    queue.memset( ptr, 0, len * sizeof(T) ) );
}



inline void sycl_device_sync( ::sycl::queue& queue ) {
  GAUXC_SYCL_ERROR( "SYCL Device Sync Failed", queue.wait_and_throw() );
}



/// Largest work-group a kernel may be launched with on the queue's device.
/// SYCL has no per-kernel query prior to bundle specialization, so this is the
/// device-wide limit clamped by the GauXC block-size convention.
inline size_t sycl_kernel_max_threads_per_block( const ::sycl::queue& queue ) {
  return queue.get_device().get_info< ::sycl::info::device::max_work_group_size >();
}

}
}

#endif
