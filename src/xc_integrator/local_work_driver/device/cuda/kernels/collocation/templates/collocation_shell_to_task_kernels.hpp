#pragma once
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include "device/common/shell_to_task.hpp"
#include <cassert>

namespace GauXC {

$py(do_grad = 'gradient' in type)\
$py(nt = 512)\

__global__ __launch_bounds__($(nt),2) void collocation_device_shell_to_task_kernel_$(type)_$(L)(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {


  __shared__ double alpha[$(nt//32)][detail::shell_nprim_max]; 
  __shared__ double coeff[$(nt//32)][detail::shell_nprim_max];
  double* my_alpha = alpha[threadIdx.x/32];
  double* my_coeff = coeff[threadIdx.x/32];

  for( auto ish = blockIdx.z; ish < nshell; ish += gridDim.z ) {
  const uint32_t ntasks      = shell_to_task[ish].ntask;
  const auto shell           = shell_to_task[ish].shell_device;
  const auto task_idx        = shell_to_task[ish].task_idx_device;
  const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;


  // Load Shell Data into registers / SM
  const uint32_t nprim = shell->nprim();
  const double3 O  = *reinterpret_cast<const double3*>(shell->O_data());

  const int global_warp_id = (threadIdx.x + blockIdx.x*blockDim.x) / cuda::warp_size;
  const int nwarp_global   = max((blockDim.x*gridDim.x) / cuda::warp_size,1);

  // Read in coeffs/exps into SM on first warp
  {
    auto* coeff_gm = shell->coeff_data();
    auto* alpha_gm = shell->alpha_data();
    static_assert( detail::shell_nprim_max == cuda::warp_size );
    const int warp_rank = threadIdx.x % cuda::warp_size;
    my_alpha[warp_rank] = alpha_gm[warp_rank];
    my_coeff[warp_rank] = coeff_gm[warp_rank];
  }

  // Loop over tasks assigned to shells
  // Place each task on a different warp + schedule across blocks
  for( int itask = global_warp_id; itask < ntasks; itask += nwarp_global ) {

    const auto*              task   = device_tasks + task_idx[itask];
    const auto* __restrict__ points_x = task->points_x;
    const auto* __restrict__ points_y = task->points_y;
    const auto* __restrict__ points_z = task->points_z;
    const uint32_t           npts   = task->npts;
    const size_t             shoff  = task_shell_offs[itask] * npts;

    auto* __restrict__ basis_eval = task->bf + shoff;
$if( do_grad )\
    auto* __restrict__ basis_x_eval = task->dbfx + shoff;
    auto* __restrict__ basis_y_eval = task->dbfy + shoff;
    auto* __restrict__ basis_z_eval = task->dbfz + shoff;
$endif\

    // Loop over points in task
    // Assign each point to separate thread within the warp
    #pragma unroll 1
    for( int ipt = threadIdx.x % cuda::warp_size; ipt < npts; ipt += cuda::warp_size ) {
      //const double3 point = points[ipt];
      double3 point;
      point.x = points_x[ipt];
      point.y = points_y[ipt];
      point.z = points_z[ipt];


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

      #pragma unroll 1
      for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = my_alpha[i];
        const auto e = my_coeff[i] * std::exp( - a * rsq );

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
        
  } // Loop over shells
} // end kernel

} // namespace GauXC
