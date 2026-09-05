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
#include "exceptions/onemkl_exception.hpp"
#include "sycl_util.hpp"

#ifdef GAUXC_HAS_SYCL
#include <memory>

namespace GauXC {
namespace util  {

/// oneMKL takes the SYCL queue directly rather than an opaque handle, so the
/// GauXC BLAS "handle" is a thin non-owning binding to the queue the BLAS
/// calls are to be issued on. This keeps the DeviceBackend API identical to
/// the cuBLAS / hipBLAS backends.
struct onemkl_handle {

  std::shared_ptr<sycl_queue> queue_ptr;

  inline onemkl_handle() = delete;
  inline onemkl_handle( std::shared_ptr<sycl_queue> q ) :
    queue_ptr(std::move(q)) { }

  inline ~onemkl_handle() noexcept = default;

  onemkl_handle( const onemkl_handle& ) = delete;
  onemkl_handle( onemkl_handle&& ) noexcept = default;

  inline ::sycl::queue& queue() { return queue_ptr->queue; }
  inline operator ::sycl::queue&() { return queue_ptr->queue; }

};


inline static ::sycl::queue& get_queue( onemkl_handle& handle ) {
  return handle.queue();
}

}
}

#endif
