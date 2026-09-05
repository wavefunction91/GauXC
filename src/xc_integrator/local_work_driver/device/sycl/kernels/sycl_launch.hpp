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
#include <sycl/sycl.hpp>
#include "device_specific/sycl_device_constants.hpp"
#include "exceptions/sycl_exception.hpp"

namespace GauXC {
namespace sycl  {

/// Launch geometry expressed the way the CUDA backend expresses it, so that
/// the ported kernels keep their original (x fastest) block/grid arithmetic.
/// SYCL orders nd_range dimensions the other way round, which the conversions
/// below take care of: CUDA .x -> dim 2, .y -> dim 1, .z -> dim 0.
struct launch_dim {
  size_t x, y, z;
  launch_dim( size_t x_ = 1, size_t y_ = 1, size_t z_ = 1 ) :
    x(x_), y(y_), z(z_) { }
};

inline ::sycl::range<3> local_range( launch_dim threads ) {
  return ::sycl::range<3>( threads.z, threads.y, threads.x );
}

inline ::sycl::range<3> global_range( launch_dim blocks, launch_dim threads ) {
  return ::sycl::range<3>( blocks.z * threads.z, blocks.y * threads.y,
                           blocks.x * threads.x );
}

inline ::sycl::nd_range<3> nd_range( launch_dim blocks, launch_dim threads ) {
  return ::sycl::nd_range<3>( global_range(blocks,threads), local_range(threads) );
}

/// CUDA-style work-item queries.
///
/// The kernels in this backend are ports of CUDA kernels that read threadIdx /
/// blockIdx / blockDim / gridDim from anywhere in the call tree, including
/// device functions several levels below the kernel entry point. Threading an
/// nd_item argument through every one of those signatures is noise, so the
/// item is recovered where it is used via the free-function queries from
/// sycl_ext_oneapi_free_function_queries.
///
/// Valid only inside a kernel launched over an nd_range<3> -- which is what
/// launch_kernel below always does.
inline ::sycl::nd_item<3> this_item() {
  return ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
}

inline ::sycl::sub_group this_subgroup() {
  return ::sycl::ext::oneapi::this_work_item::get_sub_group();
}

inline ::sycl::group<3> this_group() {
  return ::sycl::ext::oneapi::this_work_item::get_work_group<3>();
}

/// Work-group local memory, declared at the point of use the way __shared__ is
/// in the CUDA kernels. T is the full array type, e.g.
///
///   auto& sm = local_mem< double[4][warp_size][BLOCK+1] >();
///
/// The extent must be a compile-time constant, as it is for every __shared__
/// declaration in the CUDA backend. Contents are uninitialized on entry, which
/// matches __shared__; the kernels that need zeros write them explicitly.
template <typename T>
inline T& local_mem() {
  return *::sycl::ext::oneapi::group_local_memory_for_overwrite<T>( this_group() );
}

/// Enqueue a kernel over a CUDA-style launch configuration. The sub-group size
/// is pinned so that the ported sub-group reductions keep their CUDA warp
/// semantics.
///
/// The kernel op takes no arguments: it recovers its work-item coordinates via
/// this_item() and its local memory via local_mem(), exactly as a CUDA kernel
/// reads threadIdx and declares __shared__.
template <typename KernelOp>
void launch_kernel( ::sycl::queue& stream, launch_dim blocks,
  launch_dim threads, const char* name, KernelOp op ) {

  if( !blocks.x or !blocks.y or !blocks.z ) return;

  GAUXC_SYCL_ERROR( std::string(name) + " Launch Failed",
    stream.parallel_for( nd_range(blocks, threads),
      [=](::sycl::nd_item<3>) [[sycl::reqd_sub_group_size(warp_size)]] {
        op();
      })
  );

}

}
}
