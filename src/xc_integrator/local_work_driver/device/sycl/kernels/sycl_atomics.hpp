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

namespace GauXC {

/// Device-scope relaxed atomic add on global memory. Matches the ordering
/// guarantees CUDA's atomicAdd provides to the kernels ported from the CUDA
/// backend (i.e. none beyond atomicity of the update itself).
template <typename T>
inline void atomic_add_device( T* addr, T val ) {
  ::sycl::atomic_ref< T, ::sycl::memory_order::relaxed,
                      ::sycl::memory_scope::device,
                      ::sycl::access::address_space::global_space > ref( *addr );
  ref.fetch_add( val );
}

/// Work-group-scope relaxed atomic add on work-group local memory
template <typename T>
inline void atomic_add_local( T* addr, T val ) {
  ::sycl::atomic_ref< T, ::sycl::memory_order::relaxed,
                      ::sycl::memory_scope::work_group,
                      ::sycl::access::address_space::local_space > ref( *addr );
  ref.fetch_add( val );
}

/// Work-group-scope relaxed atomic post-increment on work-group local memory
template <typename T>
inline T atomic_fetch_add_local( T* addr, T val ) {
  ::sycl::atomic_ref< T, ::sycl::memory_order::relaxed,
                      ::sycl::memory_scope::work_group,
                      ::sycl::access::address_space::local_space > ref( *addr );
  return ref.fetch_add( val );
}

}
