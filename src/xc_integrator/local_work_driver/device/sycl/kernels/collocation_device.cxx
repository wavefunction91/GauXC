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
#include "device_specific/sycl_util.hpp"
#include "exceptions/sycl_exception.hpp"
#include <gauxc/xc_task.hpp>
#include <gauxc/shell.hpp>

#include "device/common/collocation_device.hpp"
#include "device/sycl/kernels/collocation_masked_kernels.hpp"
#include "device/sycl/kernels/collocation_masked_combined_kernels.hpp"
#include "device/sycl/kernels/collocation_shell_to_task_kernels.hpp"

#include "device_specific/sycl_device_constants.hpp"
#include "device/sycl/kernels/sycl_launch.hpp"

#include <algorithm>

#define GAUXC_SYCL_MAX_L 4

namespace GauXC {

namespace detail {

/// Common launcher for the shell-to-task collocation kernels. The CUDA
/// backend encodes its launch geometry in __launch_bounds__; the shell
/// scratch that used to live in a launch-time local_accessor is now
/// declared inside each kernel body via local_mem(), the way __shared__ is
/// declared in the CUDA original.
///
/// nthreads: work-group size the kernel was tuned for
/// nwarp:    rows of shell scratch, i.e. nthreads / warp_size (kept only to
///           check the launch geometry -- the kernels declare their own
///           local_mem<> extent internally)
template <uint32_t nthreads, uint32_t nwarp, typename KernelOp, typename... Args>
void launch_shell_to_task( ::sycl::queue& queue, uint32_t ntask_average,
  uint32_t nshells, KernelOp op, Args&&... args ) {

  static_assert( nthreads == nwarp * sycl::warp_size,
    "Shell scratch must have one row per sub-group" );

  if( !nshells ) return;

  const uint32_t nwarp_per_block = nthreads / sycl::warp_size;
  const uint32_t n_task_blocks   = util::div_ceil( ntask_average, nwarp_per_block );

  // dim 0 -> shells (CUDA z), dim 2 -> points/tasks (CUDA x)
  ::sycl::range<3> global( nshells, 1, n_task_blocks * nthreads );
  ::sycl::range<3> local ( 1, 1, nthreads );

  GAUXC_SYCL_ERROR( "Collocation Shell-To-Task Launch Failed",
    queue.parallel_for( ::sycl::nd_range<3>(global, local),
      [=](::sycl::nd_item<3>)
      [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
        op( nshells, args... );
      })
  );

}

}

 
template <typename T>
void eval_collocation_masked(
  size_t            nshells,
  size_t            nbf,
  size_t            npts,
  const Shell<T>*   shells_device,
  const size_t*     mask_device,
  const size_t*     offs_device,
  const T*          pts_device,
  T*                eval_device,
  device_queue queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  const auto nmax_threads = util::sycl_kernel_max_threads_per_block( stream );
  const auto max_warps_per_thread_block = nmax_threads / sycl::warp_size;

  ::sycl::range<2> local( max_warps_per_thread_block, sycl::warp_size );
  ::sycl::range<2> global(
    util::div_ceil( nshells, local[0] ) * local[0],
    util::div_ceil( npts,    local[1] ) * local[1] );

  GAUXC_SYCL_ERROR( "Collocation Masked Launch Failed",
    stream.parallel_for( ::sycl::nd_range<2>(global, local),
      [=](::sycl::nd_item<2>)
      [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
        collocation_device_masked_kernel<T>( nshells, nbf, npts,
          shells_device, mask_device, offs_device, pts_device, eval_device );
      })
  );

}
 
template             
void eval_collocation_masked(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  device_queue    queue
);




template <typename T>
void eval_collocation_masked_combined(
  size_t            ntasks,
  size_t            npts_max,
  size_t            nshells_max,
  Shell<T>*         shells_device,
  XCDeviceTask*     device_tasks,
  device_queue queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  const auto nmax_threads = util::sycl_kernel_max_threads_per_block( stream );
  const auto max_warps_per_thread_block = nmax_threads / sycl::warp_size;

  ::sycl::range<3> local( 1, max_warps_per_thread_block, sycl::warp_size );
  ::sycl::range<3> global( ntasks,
    util::div_ceil( nshells_max, local[1] ) * local[1],
    util::div_ceil( npts_max,    local[2] ) * local[2] );

  GAUXC_SYCL_ERROR( "Collocation Masked Combined Launch Failed",
    stream.parallel_for( ::sycl::nd_range<3>(global, local),
      [=](::sycl::nd_item<3>)
      [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
        collocation_device_masked_combined_kernel<T>( ntasks,
          shells_device, device_tasks );
      })
  );
     
}

template
void eval_collocation_masked_combined(
  size_t            ntasks,
  size_t            npts_max,
  size_t            nshells_max,
  Shell<double>*    shells_device,
  XCDeviceTask*     device_tasks,
  device_queue queue
);




template <typename T>
void eval_collocation_masked_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   mask_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z,
  device_queue queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  const auto nmax_threads = util::sycl_kernel_max_threads_per_block( stream );
  const auto max_warps_per_thread_block = nmax_threads / sycl::warp_size;

  ::sycl::range<2> local( max_warps_per_thread_block, sycl::warp_size );
  ::sycl::range<2> global(
    util::div_ceil( nshells, local[0] ) * local[0],
    util::div_ceil( npts,    local[1] ) * local[1] );

  GAUXC_SYCL_ERROR( "Collocation Masked Deriv1 Launch Failed",
    stream.parallel_for( ::sycl::nd_range<2>(global, local),
      [=](::sycl::nd_item<2>)
      [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
        collocation_device_masked_kernel_deriv1<T>( nshells, nbf, npts,
          shells_device, mask_device, offs_device, pts_device, eval_device,
          deval_device_x, deval_device_y, deval_device_z );
      })
  );

}

template
void eval_collocation_masked_deriv1(
  size_t               nshells,
  size_t               nbf,
  size_t               npts,
  const Shell<double>* shells_device,
  const size_t*        mask_device,
  const size_t*        offs_device,
  const double*        pts_device,
  double*              eval_device,
  double*              deval_device_x,
  double*              deval_device_y,
  double*              deval_device_z,
  device_queue    queue
);




template <typename T>
void eval_collocation_masked_combined_deriv1(
  size_t        ntasks,
  size_t        npts_max,
  size_t        nshells_max,
  Shell<T>*     shells_device,
  XCDeviceTask* device_tasks,
  device_queue queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  const auto nmax_threads = util::sycl_kernel_max_threads_per_block( stream );

  ::sycl::range<3> local( 1, nmax_threads / sycl::warp_size, sycl::warp_size );
  ::sycl::range<3> global( ntasks,
    util::div_ceil( nshells_max, local[1] ) * local[1],
    util::div_ceil( npts_max,    local[2] ) * local[2] );

  GAUXC_SYCL_ERROR( "Collocation Masked Combined Deriv1 Launch Failed",
    stream.parallel_for( ::sycl::nd_range<3>(global, local),
      [=](::sycl::nd_item<3>)
      [[sycl::reqd_sub_group_size(GauXC::sycl::warp_size)]] {
        collocation_device_masked_combined_kernel_deriv1<T>( ntasks,
          shells_device, device_tasks );
      })
  );
     
}

template
void eval_collocation_masked_combined_deriv1(
  size_t                ntasks,
  size_t                npts_max,
  size_t                nshells_max,
  Shell<double>*        shells_device,
  XCDeviceTask* device_tasks,
  device_queue queue
);


template <typename... Args>
void dispatch_shell_to_task_collocation( ::sycl::queue& queue, int32_t l,
  bool pure, uint32_t ntask_average, uint32_t nshells, Args&&... args ) {

  // Launch geometry mirrors the CUDA backend: one sub-group per task,
  // one work-group row of shell scratch per sub-group, shells over dim 0
  switch(l) {
    case 0:
      detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
        [=](uint32_t nsh, auto... a) {
          collocation_device_shell_to_task_kernel_cartesian_0( nsh, a... );
        }, std::forward<Args>(args)... );
      break;
    case 1:
      if( pure )
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_1( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_1( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 2:
      if( pure )
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_2( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_2( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 3:
      if( pure )
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_3( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_3( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 4:
      if( pure )
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_4( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_4( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    default: GAUXC_GENERIC_EXCEPTION("SYCL L_MAX = 4");
  }
}


void eval_collocation_shell_to_task(
  uint32_t                    max_l,
  AngularMomentumShellToTaskBatch* l_batched_shell_to_task,
  XCDeviceTask*               device_tasks,
  device_queue           queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  for( auto l = 0u; l <= max_l; ++l ) {
    auto pure = l_batched_shell_to_task[l].pure;
    auto shell_to_task_device = l_batched_shell_to_task[l].shell_to_task_device;
    auto nshells = l_batched_shell_to_task[l].nshells_in_batch;
    auto ntask_average = std::max(1ul, l_batched_shell_to_task[l].ntask_average);
    dispatch_shell_to_task_collocation( stream, l, pure, ntask_average,
      nshells, shell_to_task_device, device_tasks );
  }

}


template <typename... Args>
void dispatch_shell_to_task_collocation_gradient( ::sycl::queue& queue, int32_t l,
  bool pure, uint32_t ntask_average, uint32_t nshells, Args&&... args ) {

  // Launch geometry mirrors the CUDA backend: one sub-group per task,
  // one work-group row of shell scratch per sub-group, shells over dim 0
  switch(l) {
    case 0:
      detail::launch_shell_to_task<512,16>( queue, ntask_average, nshells,
        [=](uint32_t nsh, auto... a) {
          collocation_device_shell_to_task_kernel_cartesian_gradient_0( nsh, a... );
        }, std::forward<Args>(args)... );
      break;
    case 1:
      if( pure )
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_gradient_1( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_gradient_1( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 2:
      if( pure )
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_gradient_2( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_gradient_2( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 3:
      if( pure )
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_gradient_3( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_gradient_3( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 4:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_gradient_4( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_gradient_4( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    default: GAUXC_GENERIC_EXCEPTION("SYCL L_MAX = 4");
  }
}


void eval_collocation_shell_to_task_gradient(
  uint32_t                    max_l,
  AngularMomentumShellToTaskBatch* l_batched_shell_to_task,
  XCDeviceTask*               device_tasks,
  device_queue           queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  for( auto l = 0u; l <= max_l; ++l ) {
    auto pure = l_batched_shell_to_task[l].pure;
    auto shell_to_task_device = l_batched_shell_to_task[l].shell_to_task_device;
    auto nshells = l_batched_shell_to_task[l].nshells_in_batch;
    auto ntask_average = std::max(1ul, l_batched_shell_to_task[l].ntask_average);
    dispatch_shell_to_task_collocation_gradient( stream, l, pure, ntask_average,
      nshells, shell_to_task_device, device_tasks );
  }

}


template <typename... Args>
void dispatch_shell_to_task_collocation_hessian( ::sycl::queue& queue, int32_t l,
  bool pure, uint32_t ntask_average, uint32_t nshells, Args&&... args ) {

  // Launch geometry mirrors the CUDA backend: one sub-group per task,
  // one work-group row of shell scratch per sub-group, shells over dim 0
  switch(l) {
    case 0:
      detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
        [=](uint32_t nsh, auto... a) {
          collocation_device_shell_to_task_kernel_cartesian_hessian_0( nsh, a... );
        }, std::forward<Args>(args)... );
      break;
    case 1:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_hessian_1( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_hessian_1( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 2:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_hessian_2( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_hessian_2( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 3:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_hessian_3( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_hessian_3( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 4:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_hessian_4( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_hessian_4( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    default: GAUXC_GENERIC_EXCEPTION("SYCL L_MAX = 4");
  }
}


void eval_collocation_shell_to_task_hessian(
  uint32_t                    max_l,
  AngularMomentumShellToTaskBatch* l_batched_shell_to_task,
  XCDeviceTask*               device_tasks,
  device_queue           queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  for( auto l = 0u; l <= max_l; ++l ) {
    auto pure = l_batched_shell_to_task[l].pure;
    auto shell_to_task_device = l_batched_shell_to_task[l].shell_to_task_device;
    auto nshells = l_batched_shell_to_task[l].nshells_in_batch;
    auto ntask_average = std::max(1ul, l_batched_shell_to_task[l].ntask_average);
    dispatch_shell_to_task_collocation_hessian( stream, l, pure, ntask_average,
      nshells, shell_to_task_device, device_tasks );
  }

}


template <typename... Args>
void dispatch_shell_to_task_collocation_laplacian( ::sycl::queue& queue, int32_t l,
  bool pure, uint32_t ntask_average, uint32_t nshells, Args&&... args ) {

  // Launch geometry mirrors the CUDA backend: one sub-group per task,
  // one work-group row of shell scratch per sub-group, shells over dim 0
  switch(l) {
    case 0:
      detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
        [=](uint32_t nsh, auto... a) {
          collocation_device_shell_to_task_kernel_cartesian_laplacian_0( nsh, a... );
        }, std::forward<Args>(args)... );
      break;
    case 1:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_laplacian_1( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_laplacian_1( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 2:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_laplacian_2( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_laplacian_2( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 3:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_laplacian_3( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_laplacian_3( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 4:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_laplacian_4( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_laplacian_4( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    default: GAUXC_GENERIC_EXCEPTION("SYCL L_MAX = 4");
  }
}


void eval_collocation_shell_to_task_laplacian(
  uint32_t                    max_l,
  AngularMomentumShellToTaskBatch* l_batched_shell_to_task,
  XCDeviceTask*               device_tasks,
  device_queue           queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  for( auto l = 0u; l <= max_l; ++l ) {
    auto pure = l_batched_shell_to_task[l].pure;
    auto shell_to_task_device = l_batched_shell_to_task[l].shell_to_task_device;
    auto nshells = l_batched_shell_to_task[l].nshells_in_batch;
    auto ntask_average = std::max(1ul, l_batched_shell_to_task[l].ntask_average);
    dispatch_shell_to_task_collocation_laplacian( stream, l, pure, ntask_average,
      nshells, shell_to_task_device, device_tasks );
  }

}


template <typename... Args>
void dispatch_shell_to_task_collocation_lapgrad( ::sycl::queue& queue, int32_t l,
  bool pure, uint32_t ntask_average, uint32_t nshells, Args&&... args ) {

  // Launch geometry mirrors the CUDA backend: one sub-group per task,
  // one work-group row of shell scratch per sub-group, shells over dim 0
  switch(l) {
    case 0:
      detail::launch_shell_to_task<256,8>( queue, ntask_average, nshells,
        [=](uint32_t nsh, auto... a) {
          collocation_device_shell_to_task_kernel_cartesian_lapgrad_0( nsh, a... );
        }, std::forward<Args>(args)... );
      break;
    case 1:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_lapgrad_1( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_lapgrad_1( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 2:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_lapgrad_2( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_lapgrad_2( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 3:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_lapgrad_3( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_lapgrad_3( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    case 4:
      if( pure )
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_spherical_lapgrad_4( nsh, a... );
          }, std::forward<Args>(args)... );
      else
        detail::launch_shell_to_task<128,4>( queue, ntask_average, nshells,
          [=](uint32_t nsh, auto... a) {
            collocation_device_shell_to_task_kernel_cartesian_lapgrad_4( nsh, a... );
          }, std::forward<Args>(args)... );
      break;
    default: GAUXC_GENERIC_EXCEPTION("SYCL L_MAX = 4");
  }
}


void eval_collocation_shell_to_task_lapgrad(
  uint32_t                    max_l,
  AngularMomentumShellToTaskBatch* l_batched_shell_to_task,
  XCDeviceTask*               device_tasks,
  device_queue           queue
) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  for( auto l = 0u; l <= max_l; ++l ) {
    auto pure = l_batched_shell_to_task[l].pure;
    auto shell_to_task_device = l_batched_shell_to_task[l].shell_to_task_device;
    auto nshells = l_batched_shell_to_task[l].nshells_in_batch;
    auto ntask_average = std::max(1ul, l_batched_shell_to_task[l].ntask_average);
    dispatch_shell_to_task_collocation_lapgrad( stream, l, pure, ntask_average,
      nshells, shell_to_task_device, device_tasks );
  }

}


} // namespace GauXC
