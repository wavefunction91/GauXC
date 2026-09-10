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
#include <oneapi/dpl/execution>
#include <oneapi/dpl/numeric>

#include "sycl_collision_detection.hpp"
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_vector_types.hpp"
#include "xc_integrator/local_work_driver/device/sycl/kernels/sycl_launch.hpp"

namespace GauXC         {
namespace load_balancer {
namespace sycl          {

// NOTE: `sycl` here names GauXC::load_balancer::sycl (this namespace), which
// also shadows GauXC::sycl (device constants, brought in via `using namespace`
// below) and the global ::sycl namespace. Every reference to the SYCL runtime
// itself must therefore be spelled with a leading `::sycl::` in this file.
using namespace GauXC::sycl;

inline
int cube_sphere_intersect(
  const double3 lo,
  const double3 up,
  const double3 center,
  const double  rad
) {

  double dist = rad * rad;

  if( center.x < lo.x ) {
    const double r_lo = center.x - lo.x;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.x > up.x ) {
    const double r_up = center.x - up.x;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;

  if( center.y < lo.y ) {
    const double r_lo = center.y - lo.y;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.y > up.y ) {
    const double r_up = center.y - up.y;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  if( dist < 0. ) return false;


  if( center.z < lo.z ) {
    const double r_lo = center.z - lo.z;
    const double dist_lo = r_lo * r_lo;
    dist -= dist_lo;
  } else if( center.z > up.z ) {
    const double r_up = center.z - up.z;
    const double dist_up = r_up * r_up;
    dist -= dist_up;
  }

  return dist > 0.;

}


void collision_detection_gpu(
          size_t ncubes,
          size_t nspheres,
          size_t LD_bit,
    const double* low_points,
    const double* high_points,
    const double* centers,
    const double* radii,
         int32_t* collisions,
         int32_t* counts
) {
  auto it = GauXC::sycl::this_item();
  const size_t nspheres_block = (nspheres + 31) / 32;
  for (int i = it.get_global_id(2); i < ncubes; i += it.get_local_range(2) * it.get_group_range(2)) {
    counts[i] = 0;
    double3 low_point;
    double3 high_point;
    low_point.x = low_points[3*i+0];
    low_point.y = low_points[3*i+1];
    low_point.z = low_points[3*i+2];

    high_point.x = high_points[3*i+0];
    high_point.y = high_points[3*i+1];
    high_point.z = high_points[3*i+2];


    for (int j_block = 0; j_block < nspheres_block; j_block++) {
      int temp_collisions = 0;
      for (int j_inner = 0; j_inner < 32; j_inner++) {
        int j = j_block * 32 + j_inner;
        if (j < nspheres) {
          double3 center;
          double radius;
          center.x = centers[3*j+0];
          center.y = centers[3*j+1];
          center.z = centers[3*j+2];

          radius = radii[j];
          temp_collisions |= (cube_sphere_intersect(low_point, high_point, center, radius) ? 1 << (j_inner) : 0);
        }
      }
      collisions[i * LD_bit + j_block] = temp_collisions;
      counts[i] += ::sycl::popcount(temp_collisions);
    }
  }
}


static constexpr int32_t buffer_size = 8;
static constexpr int32_t element_size = 32;
static constexpr int32_t buffer_size_bits = buffer_size * element_size;

// This kernel converts the bitvector produced by the collision detection kernel above into a position list.
// For simplicity, the collision detection kernel stores its output as a bitvector. However, the `shell_list`
// of the task is a list of the qualifying indexes, so we must convert the bitvector to a position list.
//
// We take this chance to compute the nbe value from the shell sizes since the data is already being read in
void bitvector_to_position_list(
           size_t  ncubes,
           size_t  nspheres,
           size_t  LD_bit,
    const int32_t* collisions,
    const int32_t* counts,
    const  size_t* shell_size,
          int32_t* position_list,
           size_t* nbe_list
) {
  auto it = GauXC::sycl::this_item();
  auto& collisions_buffer =
    GauXC::sycl::local_mem< int32_t[warp_size][warp_size][buffer_size] >();

  // We are converting a large number of small bitvectors into position lists. For this reason, I am assigning a single thread to each bitvector
  // This avoids having to do popcounts and warp wide reductions, but hurts the memory access pattern

  // All threads in a warp must be active to do shared memory loads, so we seperate out the threadId.x
  for (int i_base = it.get_local_id(1) * it.get_local_range(2) + it.get_group(2) * it.get_local_range(2) * it.get_local_range(1);
       i_base < ncubes;
       i_base += it.get_local_range(2) * it.get_local_range(1) * it.get_group_range(2)) {
    const int i = i_base + it.get_local_id(2);
    int32_t* out = position_list;
    if (i != 0 && i < ncubes) {
      out += counts[i-1];
    }

    int current = 0;
    size_t nbe = 0;
    size_t nsphere_blocks = (nspheres + buffer_size_bits - 1) / buffer_size_bits;
    for (int j_block = 0; j_block < nsphere_blocks; j_block++) {
      // Each thread has a buffer of length BUFFER_SIZE. All the threads in the warp work to
      // load this data in a coalesced way (at least as much as possible)
      for (int buffer_loop = 0; buffer_loop < warp_size; buffer_loop += warp_size/buffer_size) {
        const int t_id_x        = it.get_local_id(2) % buffer_size;
        const int buffer_thread = it.get_local_id(2) / buffer_size;
        const int buffer_idx    = buffer_thread + buffer_loop;
        if (j_block * buffer_size_bits + t_id_x * element_size < nspheres && i_base + buffer_idx < ncubes) {
          collisions_buffer[it.get_local_id(1)][buffer_idx][t_id_x] = collisions[(i_base + buffer_idx) * LD_bit + j_block * buffer_size + t_id_x];
        }
      }

      ::sycl::group_barrier( it.get_sub_group() );
      if (i < ncubes) {  // Once the data has been loaded, we exclude the threads not corresponding to a bitvector
        // We have loaded in BUFFER_SIZE_BITS elements to be processed by each warp
        for (int j_inner = 0; j_inner < buffer_size_bits && j_block * buffer_size_bits + j_inner < nspheres; j_inner++) {
          const int j = buffer_size_bits * j_block + j_inner;
          const int j_int = j_inner / element_size;
          const int j_bit = j_inner % element_size;
          if( collisions_buffer[it.get_local_id(1)][it.get_local_id(2)][j_int] & (1 << (j_bit)) ) {
            out[current++] = j;
            nbe += shell_size[j];
          }
        }
      }
      ::sycl::group_barrier( it.get_sub_group() );
    }
    if (i < ncubes) {
      nbe_list[i] = nbe;
    }
  }
}

size_t compute_scratch( size_t ncubes, int32_t* counts_device ) {
    // Computes amount of memory that will be required to do the inclusive sum
    return 0;
}

void collision_detection( size_t        ncubes,
                          size_t        nspheres,
                          size_t        LD_bit,
                          const double* low_points_device,
                          const double* high_points_device,
                          const double* centers_device,
                          const double* radii_device,
                                size_t  temp_storage_bytes,
                                 void * temp_storage_device,
                               int32_t* collisions_device,
                               int32_t* counts_device,
                          ::sycl::queue&  stream) {

    GauXC::sycl::launch_dim threads( max_threads_per_thread_block );
    GauXC::sycl::launch_dim blocks( util::div_ceil( ncubes, threads.x ) );

    GauXC::sycl::launch_kernel( stream, blocks, threads, "collision_detection_gpu",
      [=](){
        collision_detection_gpu(
          ncubes, nspheres, LD_bit,
          low_points_device, high_points_device, centers_device, radii_device,
          collisions_device, counts_device
        );
      });

    // Run inclusive prefix sum
    oneapi::dpl::inclusive_scan(
      oneapi::dpl::execution::make_device_policy(stream),
      counts_device, counts_device + ncubes, counts_device );

}

void compute_position_list(size_t         ncubes,
                           size_t         nspheres,
                           size_t         LD_bit,
                           const size_t*  shell_sizes_device,
                           const int32_t* collisions_device,
                           const int32_t* counts_device,
                                 int32_t* position_list_device,
                                  size_t* nbe_list_device,
                            ::sycl::queue&  stream) {
    GauXC::sycl::launch_dim threads( warp_size, warp_size );
    GauXC::sycl::launch_dim blocks( util::div_ceil( ncubes, threads.x * threads.y ) );

    // convert from bitvector to position list
    GauXC::sycl::launch_kernel( stream, blocks, threads, "bitvector_to_position_list",
      [=](){
        bitvector_to_position_list(
          ncubes, nspheres, LD_bit,
          collisions_device, counts_device, shell_sizes_device,
          position_list_device, nbe_list_device
        );
      });
}

}
}
}
