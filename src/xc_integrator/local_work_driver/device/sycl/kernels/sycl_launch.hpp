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

inline ::sycl::nd_item<3> this_item() {
  return ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
}

inline ::sycl::sub_group this_subgroup() {
  return ::sycl::ext::oneapi::this_work_item::get_sub_group();
}

inline ::sycl::group<3> this_group() {
  return ::sycl::ext::oneapi::this_work_item::get_work_group<3>();
}

template <typename T>
inline T& local_mem() {
  return *::sycl::ext::oneapi::group_local_memory_for_overwrite<T>( this_group() );
}

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
