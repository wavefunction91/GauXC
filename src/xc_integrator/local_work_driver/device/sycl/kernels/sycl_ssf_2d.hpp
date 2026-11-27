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
#include <numeric>

#include "sycl_extensions.hpp"
#include "sycl_atomics.hpp"
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_vector_types.hpp"
#include "common/integrator_constants.hpp"

inline constexpr static auto eps_d = std::numeric_limits<double>::epsilon();

namespace GauXC {

inline double gFrisch(double x) {
  // Frisch partition functions
  const double s_x  = x * 1.5625;
  const double s_x2 = s_x  * s_x;
  const double s_x3 = s_x  * s_x2;
  const double s_x5 = s_x3 * s_x2;
  const double s_x7 = s_x5 * s_x2;

  return ((35.) *(s_x - s_x3) + (21.) *s_x5 - (5.) *s_x7);
}


inline double sFrisch(double x) {

    if( ::sycl::fabs(x) < integrator::magic_ssf_factor<> ) return (0.5 - (0.5/ 16.0) * gFrisch(x));
    else if( x >= integrator::magic_ssf_factor<> ) return 0.;
    else                               return 1.;
}

template <uint32_t weight_unroll,             // Unrolling factor
          uint32_t weight_thread_block,       // Number of work-items / work-group
          uint32_t weight_thread_block_per_sm // Work-groups / XeCore
          >
void modify_weights_ssf_kernel_2d(
  int32_t npts, int32_t natoms,
  const double* RAB, int32_t ldRAB, const double* coords, const double* dist,
  size_t lddist, const int32_t* iparent, const double* dist_nearest,
  double* weights ) {

  static_assert( weight_unroll == 4 );

  auto it = ::sycl::ext::oneapi::this_work_item::get_nd_item<2>();
  constexpr uint32_t warps_per_block = weight_thread_block / sycl::warp_size;
  auto& jCounter_sm = *::sycl::ext::oneapi::group_local_memory_for_overwrite<int[warps_per_block]>(
    ::sycl::ext::oneapi::this_work_item::get_work_group<2>() );

  constexpr double weight_tol = integrator::ssf_weight_tol;

  const auto sg = it.get_sub_group();

  const auto local_x = it.get_local_id(1);
  const auto local_y = it.get_local_id(0);
  const auto range_x = it.get_local_range(1);

  int natom_block = ((natoms + range_x - 1) / range_x) * range_x;

  const size_t tid_x = it.get_global_id(0);
  const size_t nt_x  = it.get_global_range(0);

  int* jCounter = &jCounter_sm[local_y];

  // Each sub-group will work together on a point
  for( size_t ipt = tid_x; ipt < (size_t)npts; ipt += nt_x ) {

    const auto iParent = iparent[ipt];

    double sum = 0.; 
    double parent_weight = 0.;

    const double* const local_dist_scratch = dist + ipt * lddist;
    const double dist_cutoff = 0.5 * (1 - integrator::magic_ssf_factor<> ) * 
      dist_nearest[ipt];
    if( local_dist_scratch[iParent] < dist_cutoff ) continue;

    // Do iParent First
    {

      const double ri = local_dist_scratch[ iParent ];
      const double* const local_rab = RAB + iParent * ldRAB;

      parent_weight = 1.;
      for( int jCenter = local_x; jCenter < natom_block; jCenter += range_x ) {
        double contribution = 1.0;
        if (jCenter < natoms && iParent != jCenter) {
          const double rj = local_dist_scratch[ jCenter ];
          const double mu = (ri - rj) * local_rab[ jCenter ]; // XXX: RAB is symmetric
          contribution = sFrisch( mu );
        }
        contribution = sycl::warp_reduce_prod<sycl::warp_size>(contribution);
        contribution = ::sycl::group_broadcast(sg, contribution, 0);

        parent_weight *= contribution;

        if (parent_weight < weight_tol) break;
      }
    }

    if( parent_weight < eps_d ) {
      if (local_x == 0)
        weights[ipt] = 0.;
      ::sycl::group_barrier(sg);
      continue;
    }

    // Initialize each counter to 0
    if (local_x == 0) {
      jCounter[0] = 0;
    }
    ::sycl::group_barrier(sg);

    // Each work-item will process an iCenter. Atomic operations are used to
    // assign an iCenter value to each work-item.
    int iCenter = atomic_fetch_add_local(jCounter, 1);
    if (iCenter >= iParent) iCenter++; // iCenter == iParent is skipped

    // The entire sub-group processes the same jCenter value at the same time
    int jCenter = 0;

    const double* local_rab = RAB + iCenter * ldRAB;
    double ri = local_dist_scratch[ iCenter ];
    double ps = 1.;
    int iCount = 0; 
    bool cont = (iCenter < natoms);

    // We will continue iterating until all of the work-items have cont == false
    while (::sycl::any_of_group(sg, cont)) {
      if (cont) {
        double2 rj[weight_unroll/2];
        double2 rab_val[weight_unroll/2];
        double mu[weight_unroll];
        iCount += weight_unroll;

        #pragma unroll
        for (int k = 0; k < weight_unroll/2; k++) {
          rj[k]      = *((const double2*)(local_dist_scratch + jCenter) + k);
          rab_val[k] = *((const double2*)(local_rab          + jCenter) + k); 
        }

        #pragma unroll
        for (int k = 0; k < weight_unroll/2; k++) {
          mu[2*k+0] = (ri - rj[k].x) * rab_val[k].x; // XXX: RAB is symmetric
          mu[2*k+1] = (ri - rj[k].y) * rab_val[k].y; 
        }

        #pragma unroll
        for (int k = 0; k < weight_unroll; k++) {
          if((iCenter != jCenter + k) && (jCenter + k < natoms)) {
            mu[k] = sFrisch( mu[k] );
            ps *= mu[k];
          }
        }

        // A work-item is done with an iCenter based on 2 conditions: weight
        // tolerance, or if it has seen all of the jCenters
        if( !(ps > weight_tol && iCount < (int)lddist )) {
          // In that case the work-item begins processing another iCenter
          sum += ps;
          iCenter = atomic_fetch_add_local(jCounter, 1);
          if (iCenter >= iParent) iCenter++;

          // If there are no more iCenters left to process, signal ready to exit
          cont = (iCenter < natoms);
          ri = local_dist_scratch[ iCenter ];
          local_rab = RAB + iCenter * ldRAB;
          ps = 1.;
          iCount = 0;
        }
      }
      // Wraps jCenter around. This was faster than modulo
      jCenter += weight_unroll;
      jCenter = (jCenter < ldRAB) ? jCenter : 0;
    }

    // All of the work-items then sum their contributions. Only lane 0 needs to
    // add the parent contribution.
    ::sycl::group_barrier(sg);
    sum = sycl::warp_reduce_sum<sycl::warp_size>(sum);
    if (local_x == 0) {
      sum += parent_weight;
      weights[ipt] *= parent_weight / sum;
    }

    ::sycl::group_barrier(sg);
  }

}


void partition_weights_ssf_2d( int32_t npts, int32_t natoms, const double* RAB,
  int32_t ldRAB, const double* coords, const double* dist, int32_t lddist,
  const int32_t* iparent, const double* dist_nearest, double* weights,
  ::sycl::queue& stream );


}
