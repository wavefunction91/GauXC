#pragma once
#include <gauxc/shell.hpp>
#include <cooperative_groups.h>
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include <cassert>

namespace GauXC {

$py(do_grad = 'gradient' in type)\
__global__ __launch_bounds__(1024,1) void collocation_device_shell_to_task_kernel_$(type)_$(L)(
  size_t                            ntasks,
  const Shell<double>* __restrict__ shell,
  const int32_t*       __restrict__ task_idx,
  const int32_t*       __restrict__ task_shell_offs,
  XCDeviceTask*        __restrict__ device_tasks
) {

  // Load Shell Data into registers / SM
  const auto nprim = shell->nprim();
  const double3 O  = *reinterpret_cast<const double3*>(shell->O_data());

#if 0
  namespace cg = cooperative_groups; 
  auto launch_grid  = cg::this_grid();
  auto thread_block = cg::this_thread_block();
  auto block_warp_partition = cg::tiled_partition<cuda::warp_size>(thread_block);

  const auto warp_rank       = block_warp_partition.thread_rank();
  const auto block_warp_id   = block_warp_partition.meta_group_rank();
  const auto global_warp_id  = launch_grid.thread_rank() / cuda::warp_size;
  const auto nwarp_global    = launch_grid.size() / cuda::warp_size;
#else
  const int warp_rank      = threadIdx.x % 32;
  const int block_warp_id  = threadIdx.x / 32;
  const int global_warp_id = (threadIdx.x + blockIdx.x*blockDim.x) / 32;
  const int nwarp_global   = (blockDim.x*gridDim.x) / 32;
#endif

  __shared__ double alpha[detail::shell_nprim_max], coeff[detail::shell_nprim_max];
  // Read in coeffs/exps into SM on first warp
  if( block_warp_id == 0 ) {
    auto* coeff_gm = shell->coeff_data();
    auto* alpha_gm = shell->alpha_data();
    for( int i = warp_rank; i < detail::shell_nprim_max; i += cuda::warp_size ) {
       alpha[i] = alpha_gm[i];
       coeff[i] = coeff_gm[i];
    }
  }
  #if 0
  thread_block.sync(); // Sync once SM is populated 
  #else
  __syncthreads();
  #endif

  // Loop over tasks assigned to shells
  // Place each task on a different warp + schedule across blocks
  for( int itask = global_warp_id; itask < ntasks; itask += nwarp_global ) {

    const auto*              task   = device_tasks + task_idx[itask];
    const auto* __restrict__ points = reinterpret_cast<double3*>(task->points);
    const auto               npts   = task->npts;
    const auto               shoff  = task_shell_offs[itask];

    auto* __restrict__ basis_eval = task->bf + shoff*npts;
$if( do_grad )\
    auto* __restrict__ basis_x_eval = task->dbfx + shoff*npts;
    auto* __restrict__ basis_y_eval = task->dbfy + shoff*npts;
    auto* __restrict__ basis_z_eval = task->dbfz + shoff*npts;
$endif\

    // Loop over points in task
    // Assign each point to separate thread within the warp
    for( int ipt = warp_rank; ipt < npts; ipt += cuda::warp_size ) {
      const double3 point = points[ipt];

      const auto x = point.x - O.x;
      const auto y = point.y - O.y;
      const auto z = point.z - O.z;
      const auto rsq = x*x + y*y + z*z;

      // Evaluate radial part of bfn
      double radial_eval = 0.;
$if( do_grad )\
      double radial_eval_x = 0.;
      double radial_eval_y = 0.;
      double radial_eval_z = 0.;
$endif\

      for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = alpha[i];
        const auto e = coeff[i] * std::exp( - a * rsq );

        radial_eval += e;
$if( do_grad )\

        const auto ae = 2. * a * e;
        radial_eval_x -= ae * x;
        radial_eval_y -= ae * y;
        radial_eval_z -= ae * z;
$endif\
      }

      
      // Evaluate the angular part of bfn
      double ang_eval;
$if( do_grad )\
      double dang_eval_x, dang_eval_y, dang_eval_z;
$endif\

$for( i in range(len(eval_lines)) )\
      $(eval_lines[i])
$if( do_grad )\
      $(eval_lines_dx[i])
      $(eval_lines_dy[i])
      $(eval_lines_dz[i])

$endif\
      basis_eval  [ipt + $(i)*npts] = ang_eval;
$if( do_grad )\
      basis_x_eval[ipt + $(i)*npts] = dang_eval_x;
      basis_y_eval[ipt + $(i)*npts] = dang_eval_y;
      basis_z_eval[ipt + $(i)*npts] = dang_eval_z;
$endif\


$endfor\
      
    }

  }
}

}