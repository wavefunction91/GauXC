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
#include "device/common/symmetrize_mat.hpp"
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_util.hpp"
#include "sycl_launch.hpp"
#include <array>
#include <utility>

namespace GauXC {

void symmetrize_matrix_device( size_t N, double* A, size_t LDA ) {

  auto it = GauXC::sycl::this_item();

  constexpr uint32_t block_size = sycl::warp_size;

  auto& buffer = GauXC::sycl::local_mem<double[block_size][block_size+1]>();

  const size_t num_blocks = ((N + block_size - 1) / block_size);

  for (int i = it.get_group(2); i < num_blocks; i += it.get_group_range(2)) {
    // TODO This could be load balanced if need be
    const int i_coord = i * block_size;
    for (int j = i; j < num_blocks; j++) {
      const int j_coord = j * block_size;

      // Read in block to buffer
      // TODO These could be vector reads/writes if this becomes significant
      if (i_coord + it.get_local_id(1) < N && j_coord + it.get_local_id(2) < N) {
        buffer[it.get_local_id(1)][it.get_local_id(2)] =
          A[(i_coord + it.get_local_id(1)) * LDA + j_coord + it.get_local_id(2)];
      }
      ::sycl::group_barrier( it.get_group() );

      // Write buffer
      if (j_coord + it.get_local_id(1) < N && i_coord + it.get_local_id(2) < N) {
        if ((j_coord != i_coord || it.get_local_id(2) < it.get_local_id(1))) { // handles the diagonal block
          A[(j_coord + it.get_local_id(1)) * LDA + i_coord + it.get_local_id(2)] =
            buffer[it.get_local_id(2)][it.get_local_id(1)];
        }
      }
      ::sycl::group_barrier( it.get_group() );
    }
  }
}

void symmetrize_matrix_inc_device( size_t N, double* A, size_t LDA ) {

  auto it = GauXC::sycl::this_item();

  constexpr uint32_t block_size = sycl::warp_size;

  auto& buffer_0 = GauXC::sycl::local_mem<double[block_size][block_size+1]>();
  auto& buffer_1 = GauXC::sycl::local_mem<double[block_size][block_size+1]>();

  const size_t num_blocks = ((N + block_size - 1) / block_size);

  for (int i = it.get_group(2); i < num_blocks; i += it.get_group_range(2)) {
    // TODO This could be load balanced if need be
    const int i_coord = i * block_size;
    for (int j = i; j < num_blocks; j++) {
      const int j_coord = j * block_size;

      // Read in block to buffer
      // TODO These could be vector reads/writes if this becomes significant
      if (i_coord + it.get_local_id(1) < N && j_coord + it.get_local_id(2) < N) {
        buffer_0[it.get_local_id(1)][it.get_local_id(2)] =
          A[(i_coord + it.get_local_id(1)) * LDA + j_coord + it.get_local_id(2)];
      }
      if (j_coord + it.get_local_id(1) < N && i_coord + it.get_local_id(2) < N) {
        buffer_1[it.get_local_id(1)][it.get_local_id(2)] =
          A[(j_coord + it.get_local_id(1)) * LDA + i_coord + it.get_local_id(2)];
      }
      ::sycl::group_barrier( it.get_group() );

      buffer_0[it.get_local_id(1)][it.get_local_id(2)] +=
        buffer_1[it.get_local_id(2)][it.get_local_id(1)];
      buffer_0[it.get_local_id(1)][it.get_local_id(2)] *= 0.5;
      ::sycl::group_barrier( it.get_group() );

      // Write buffer
      if (j_coord + it.get_local_id(1) < N && i_coord + it.get_local_id(2) < N) {
        //if ((j_coord != i_coord || it.get_local_id(2) < it.get_local_id(1))) { // handles the diagonal block
          A[(j_coord + it.get_local_id(1)) * LDA + i_coord + it.get_local_id(2)] =
            buffer_0[it.get_local_id(2)][it.get_local_id(1)];
        //}
      }
      if (i_coord + it.get_local_id(1) < N && j_coord + it.get_local_id(2) < N) {
        //if ((j_coord != i_coord || it.get_local_id(2) > it.get_local_id(1))) { // handles the diagonal block
          A[(i_coord + it.get_local_id(1)) * LDA + j_coord + it.get_local_id(2)] =
            buffer_0[it.get_local_id(1)][it.get_local_id(2)];
        //}
      }
      ::sycl::group_barrier( it.get_group() );
    }
  }
}



void symmetrize_matrix( int32_t N, double* A, size_t LDA, device_queue queue ) {
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();
  const size_t num_blocks = ((N + sycl::warp_size - 1) / sycl::warp_size);
  // Warp size must equal max_warps_per_thread_block must equal 32
  sycl::launch_dim threads(sycl::warp_size, sycl::max_warps_per_thread_block, 1), blocks(num_blocks);

  sycl::launch_kernel( stream, blocks, threads, "symmetrize_matrix_device",
    [=](){
      symmetrize_matrix_device( N, A, LDA );
    });
}

void symmetrize_matrix_inc( int32_t N, double* A, size_t LDA, device_queue queue ) {
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();
  const size_t num_blocks = ((N + sycl::warp_size - 1) / sycl::warp_size);
  // Warp size must equal max_warps_per_thread_block must equal 32
  sycl::launch_dim threads(sycl::warp_size, sycl::max_warps_per_thread_block, 1), blocks(num_blocks);

  sycl::launch_kernel( stream, blocks, threads, "symmetrize_matrix_inc_device",
    [=](){
      symmetrize_matrix_inc_device( N, A, LDA );
    });
}
}
