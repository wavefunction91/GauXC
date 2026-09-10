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
#include "sycl_ssf_2d.hpp"
#include "device/sycl/sycl_aos_scheme1.hpp"
#include "device_specific/sycl_util.hpp"

namespace GauXC {

void partition_weights_ssf_2d( int32_t npts, int32_t natoms, const double* RAB,
  int32_t ldRAB, const double* coords, const double* dist, int32_t lddist,
  const int32_t* iparent, const double* dist_nearest, double* weights,
  ::sycl::queue& stream ) {

  constexpr auto weight_unroll =
    alg_constants::SyclAoSScheme1::weight_unroll;
  constexpr auto weight_thread_block =
    alg_constants::SyclAoSScheme1::weight_thread_block;
  constexpr auto weight_thread_block_per_sm =
    alg_constants::SyclAoSScheme1::weight_thread_block_per_sm;

  const auto num_sm =
    stream.get_device().get_info< ::sycl::info::device::max_compute_units >();

  constexpr auto warps_per_block = weight_thread_block / sycl::warp_size;

  ::sycl::range<2> local ( warps_per_block, sycl::warp_size );
  ::sycl::range<2> global( num_sm * weight_thread_block_per_sm * warps_per_block,
                           sycl::warp_size );

  GAUXC_SYCL_ERROR( "SSF 2D Weights Launch Failed",
    stream.submit([&](::sycl::handler& cgh) {
      cgh.parallel_for( ::sycl::nd_range<2>(global, local),
        [=](::sycl::nd_item<2>)
        [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
          modify_weights_ssf_kernel_2d< weight_unroll, weight_thread_block,
            weight_thread_block_per_sm >( npts, natoms, RAB,
              ldRAB, coords, dist, lddist, iparent, dist_nearest, weights );
        });
    })
  );

}

}
