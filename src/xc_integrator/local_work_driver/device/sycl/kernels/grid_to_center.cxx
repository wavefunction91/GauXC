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
#include <gauxc/util/div_ceil.hpp>
#include "grid_to_center.hpp"
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_util.hpp"
#include "device_specific/sycl_vector_types.hpp"

namespace GauXC {

void compute_grid_to_center_dist( int32_t npts, int32_t natoms,
  const double* coords, const double* points_x, const double* points_y, 
  const double* points_z, double* dist, int32_t lddist, ::sycl::queue& stream ) {

    const uint32_t distance_thread_y = sycl::max_warps_per_thread_block / 2;

    // dim 1 <-> CUDA y, dim 2 <-> CUDA x
    ::sycl::range<2> local( distance_thread_y, sycl::warp_size );
    ::sycl::range<2> global(
      util::div_ceil( npts, local[0] * distance_thread_y ) * local[0],
      util::div_ceil( natoms, local[1] ) * local[1] );

    // This kernel is launched over a 2D nd_range (matching the CUDA original's
    // dim3 threads(warp_size, distance_thread_y)), so it does not go through
    // the 3D-only launch_kernel/this_item()/local_mem() contract in
    // sycl_launch.hpp. It still uses the same free-function work-item query
    // style, just instantiated for 2 dimensions, and group_local_memory_for_
    // overwrite for the compile-time-sized point_buffer (mirrors the CUDA
    // __shared__ double3 point_buffer[warp_size]).
    GAUXC_SYCL_ERROR( "Grid-To-Center Launch Failed",
      stream.submit([&](::sycl::handler& cgh) {
        cgh.parallel_for( ::sycl::nd_range<2>(global, local),
          [=](::sycl::nd_item<2>)
          [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {

          auto it = ::sycl::ext::oneapi::this_work_item::get_nd_item<2>();
          auto& point_buffer = *::sycl::ext::oneapi::group_local_memory_for_overwrite<double3[sycl::warp_size]>(
            ::sycl::ext::oneapi::this_work_item::get_work_group<2>() );

          const auto tid_x = it.get_local_id(1);
          const auto tid_y = it.get_local_id(0);

          double3 coord_reg{0.,0.,0.};

          const int natoms_block = (natoms + sycl::warp_size-1) / sycl::warp_size;
          const int coords_block = (npts   + sycl::warp_size-1) / sycl::warp_size;

          const double3* coords_vec = (const double3*) coords;

          for (size_t j = it.get_group(1); j < (size_t)natoms_block; j += it.get_group_range(1)) {
            const int iAtom = j * sycl::warp_size + tid_x;
            // Load blocks into registers/local memory
            if (iAtom < natoms) {
              coord_reg = coords_vec[iAtom];
            }
            for (size_t i = it.get_group(0); i < (size_t)coords_block; i += it.get_group_range(0)) {
              const int iPt_load = i * sycl::warp_size + tid_x;
              if (iPt_load < npts) {
                point_buffer[tid_x].x = points_x[iPt_load];
                point_buffer[tid_x].y = points_y[iPt_load];
                point_buffer[tid_x].z = points_z[iPt_load];
              }
              ::sycl::group_barrier( it.get_group() );

              // do the computation
              for (uint32_t k = tid_y; k < sycl::warp_size; k += sycl::warp_size/2) {
                const int iPt_sm = k;
                const int iPt = i * sycl::warp_size + iPt_sm;
                const double rx = point_buffer[iPt_sm].x - coord_reg.x;
                const double ry = point_buffer[iPt_sm].y - coord_reg.y;
                const double rz = point_buffer[iPt_sm].z - coord_reg.z;

                if (iAtom < natoms and iPt < npts) {
                  dist[ iAtom + iPt * lddist ] = ::sycl::sqrt( rx*rx + ry*ry + rz*rz );
                }
              }
              ::sycl::group_barrier( it.get_group() );
            }
          }
        });
      })
    );

}

}
