#pragma once
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include "device/common/shell_to_task.hpp"
#include <cassert>

namespace GauXC {


__global__ __launch_bounds__(512,1) void collocation_device_shell_to_task_kernel_spherical_2(
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

      for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = alpha[i];
        const auto e = coeff[i] * std::exp( - a * rsq );

        radial_eval += e;
      }

      
      // Evaluate the angular part of bfn



      double ang_eval_0;
      double ang_eval_1;
      double ang_eval_2;
      double ang_eval_3;


      ang_eval_0 = sqrt_3*radial_eval*x*y;
      ang_eval_1 = sqrt_3*radial_eval*y*z;
      ang_eval_2 = radial_eval*(-x*x - y*y + 2*z*z)/2;
      ang_eval_3 = sqrt_3*radial_eval*x*z;
      basis_eval[ipt + 0*npts] = ang_eval_0;
      basis_eval[ipt + 1*npts] = ang_eval_1;
      basis_eval[ipt + 2*npts] = ang_eval_2;
      basis_eval[ipt + 3*npts] = ang_eval_3;

      ang_eval_0 = sqrt_3*radial_eval*(x*x - y*y)/2;
      basis_eval[ipt + 4*npts] = ang_eval_0;


    } // Loop over points within task
  } // Loop over tasks
  #endif
  } // Loop over shells
} // end kernel

} // namespace GauXC
