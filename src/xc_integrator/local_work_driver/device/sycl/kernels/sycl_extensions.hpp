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
#include "sycl_launch.hpp"

namespace GauXC {
namespace sycl  {

template <size_t warp_sz, typename T>
inline T warp_reduce_sum( T val ) {

  auto sg = this_subgroup();
  for(int i=(warp_sz/2); i>=1; i/=2)
    val += ::sycl::permute_group_by_xor( sg, val, i );

  return val;
}

template <size_t warp_sz, typename T>
inline T warp_reduce_prod( T val ) {

  auto sg = this_subgroup();
  for(int i=(warp_sz/2); i>=1; i/=2)
    val *= ::sycl::permute_group_by_xor( sg, val, i );

  return val;
}

template <size_t warp_sz, typename T>
inline T warp_reduce_max( T val ) {

  auto sg = this_subgroup();
  for(int i=(warp_sz/2); i>=1; i/=2)
    val = ::sycl::fmax( val, ::sycl::permute_group_by_xor( sg, val, i ) );

  return val;
}

}
}
