#pragma once
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include "device/common/shell_to_task.hpp"
#include <cassert>

namespace GauXC {

$py(do_grad = 'gradient' in type)\

__global__ __launch_bounds__(256,1) void collocation_device_shell_to_task_kernel_$(type)_$(L)(
  int32_t                         nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {


  __shared__ double alpha[detail::shell_nprim_max], coeff[detail::shell_nprim_max];

  for( auto ish = blockIdx.y; ish < nshell; ish += gridDim.y ) {
  const auto ntasks   = shell_to_task[ish].ntask;
  const auto shell    = shell_to_task[ish].shell_device;
  const auto task_idx = shell_to_task[ish].task_idx_device;
  const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;


  // Load Shell Data into registers / SM
  const auto nprim = shell->nprim();
  const double3 O  = *reinterpret_cast<const double3*>(shell->O_data());

  const int warp_rank      = threadIdx.x % cuda::warp_size;
  const int block_warp_id  = threadIdx.x / cuda::warp_size;
  const int global_warp_id = (threadIdx.x + blockIdx.x*blockDim.x) / cuda::warp_size;
  const int nwarp_global   = max((blockDim.x*gridDim.x) / cuda::warp_size,1);

  if(ish) __syncthreads(); // Sync to avoid invalidation of cache for warps working on a different shell
  // Read in coeffs/exps into SM on first warp
  if( block_warp_id == 0 ) {
    auto* coeff_gm = shell->coeff_data();
    auto* alpha_gm = shell->alpha_data();
    for( int i = warp_rank; i < detail::shell_nprim_max; i += cuda::warp_size ) {
       alpha[i] = alpha_gm[i];
       coeff[i] = coeff_gm[i];
    }
  }
  __syncthreads(); // Sync once SM is populated 

#if 1
  // Loop over tasks assigned to shells
  // Place each task on a different warp + schedule across blocks
  for( int itask = global_warp_id; itask < ntasks; itask += nwarp_global ) {

    const auto*              task   = device_tasks + task_idx[itask];
    const auto* __restrict__ points = reinterpret_cast<double3*>(task->points);
    const auto               npts   = task->npts;
    const auto               shoff  = task_shell_offs[itask] * npts;

    auto* __restrict__ basis_eval = task->bf + shoff;
$if( do_grad )\
    auto* __restrict__ basis_x_eval = task->dbfx + shoff;
    auto* __restrict__ basis_y_eval = task->dbfy + shoff;
    auto* __restrict__ basis_z_eval = task->dbfz + shoff;
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

$py(unroll_max = min(len(eval_lines),4))

$for( i in range(unroll_max) )\
      double ang_eval_$(i);
$endfor\

$py(unroll_loop_ceil = len(eval_lines)//unroll_max)
$py(idx_st = unroll_loop_ceil*unroll_max)\
$for( i in range(unroll_loop_ceil) )\
$for( j in range(unroll_max) )\
      ang_eval_$(j) = $(eval_lines[i*unroll_max + j]);
$endfor\
$for( j in range(unroll_max) )\
      basis_eval[ipt + $(i*unroll_max + j)*npts] = ang_eval_$(j);
$endfor\

$endfor\
$if( len(eval_lines)%unroll_max )\
$for( j in range(len(eval_lines)%unroll_max) )\
      ang_eval_$(j) = $(eval_lines[idx_st + j]);
$endfor\
$for( j in range(len(eval_lines)%unroll_max) )\
      basis_eval[ipt + $(idx_st + j)*npts] = ang_eval_$(j);
$endfor\

$endif\

$if(do_grad)\
$for( i in range(unroll_max) )\
      double dang_eval_x_$(i), dang_eval_y_$(i), dang_eval_z_$(i);
$endfor\

$for( i in range(unroll_loop_ceil) )\
$for( j in range(unroll_max) )\
      dang_eval_x_$(j) = $(eval_lines_dx[i*unroll_max + j]);
      dang_eval_y_$(j) = $(eval_lines_dy[i*unroll_max + j]);
      dang_eval_z_$(j) = $(eval_lines_dz[i*unroll_max + j]);
$endfor\
$for( j in range(unroll_max) )\
      basis_x_eval[ipt + $(i*unroll_max + j)*npts] = dang_eval_x_$(j);
      basis_y_eval[ipt + $(i*unroll_max + j)*npts] = dang_eval_y_$(j);
      basis_z_eval[ipt + $(i*unroll_max + j)*npts] = dang_eval_z_$(j);
$endfor\

$endfor\
$if( len(eval_lines)%unroll_max )\
$for( j in range(len(eval_lines)%unroll_max) )\
      dang_eval_x_$(j) = $(eval_lines_dx[idx_st + j]);
      dang_eval_y_$(j) = $(eval_lines_dy[idx_st + j]);
      dang_eval_z_$(j) = $(eval_lines_dz[idx_st + j]);
$endfor\
$for( j in range(len(eval_lines)%unroll_max) )\
      basis_x_eval[ipt + $(idx_st + j)*npts] = dang_eval_x_$(j);
      basis_y_eval[ipt + $(idx_st + j)*npts] = dang_eval_y_$(j);
      basis_z_eval[ipt + $(idx_st + j)*npts] = dang_eval_z_$(j);
$endfor\

$endif\
$endif\
    } // Loop over points within task
  } // Loop over tasks
  #endif
  } // Loop over shells
} // end kernel

} // namespace GauXC
