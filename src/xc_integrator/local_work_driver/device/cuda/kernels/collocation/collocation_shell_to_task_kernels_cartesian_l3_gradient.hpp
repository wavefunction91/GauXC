/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include "device/common/shell_to_task.hpp"
#include <cassert>

namespace GauXC {


__global__ __launch_bounds__(512,2) void collocation_device_shell_to_task_kernel_cartesian_gradient_3(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {


  __shared__ double alpha[16][detail::shell_nprim_max + 1]; 
  __shared__ double coeff[16][detail::shell_nprim_max + 1];
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
    auto* __restrict__ basis_x_eval = task->dbfx + shoff;
    auto* __restrict__ basis_y_eval = task->dbfy + shoff;
    auto* __restrict__ basis_z_eval = task->dbfz + shoff;


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
      double radial_eval_alpha = 0.;

      #pragma unroll 1
      for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = my_alpha[i];
        const auto e = my_coeff[i] * std::exp( - a * rsq );

        radial_eval += e;
        radial_eval_alpha += a * e;
      }

      radial_eval_alpha *= -2;

      

      // Evaluate basis function
      basis_eval[ipt + 0*npts] = radial_eval*x*x*x;
      basis_eval[ipt + 1*npts] = radial_eval*x*x*y;
      basis_eval[ipt + 2*npts] = radial_eval*x*x*z;
      basis_eval[ipt + 3*npts] = radial_eval*x*y*y;
      basis_eval[ipt + 4*npts] = radial_eval*x*y*z;
      basis_eval[ipt + 5*npts] = radial_eval*x*z*z;
      basis_eval[ipt + 6*npts] = radial_eval*y*y*y;
      basis_eval[ipt + 7*npts] = radial_eval*y*y*z;
      basis_eval[ipt + 8*npts] = radial_eval*y*z*z;
      basis_eval[ipt + 9*npts] = radial_eval*z*z*z;


    
      // Evaluate first derivative of bfn wrt x
      basis_x_eval[ipt + 0*npts] = x*x*(3*radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 1*npts] = x*y*(2*radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 2*npts] = x*z*(2*radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 3*npts] = y*y*(radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 4*npts] = y*z*(radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 5*npts] = z*z*(radial_eval + radial_eval_alpha*x*x);
      basis_x_eval[ipt + 6*npts] = radial_eval_alpha*x*y*y*y;
      basis_x_eval[ipt + 7*npts] = radial_eval_alpha*x*y*y*z;
      basis_x_eval[ipt + 8*npts] = radial_eval_alpha*x*y*z*z;
      basis_x_eval[ipt + 9*npts] = radial_eval_alpha*x*z*z*z;

      // Evaluate first derivative of bfn wrt y
      basis_y_eval[ipt + 0*npts] = radial_eval_alpha*x*x*x*y;
      basis_y_eval[ipt + 1*npts] = x*x*(radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 2*npts] = radial_eval_alpha*x*x*y*z;
      basis_y_eval[ipt + 3*npts] = x*y*(2*radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 4*npts] = x*z*(radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 5*npts] = radial_eval_alpha*x*y*z*z;
      basis_y_eval[ipt + 6*npts] = y*y*(3*radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 7*npts] = y*z*(2*radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 8*npts] = z*z*(radial_eval + radial_eval_alpha*y*y);
      basis_y_eval[ipt + 9*npts] = radial_eval_alpha*y*z*z*z;

      // Evaluate first derivative of bfn wrt z
      basis_z_eval[ipt + 0*npts] = radial_eval_alpha*x*x*x*z;
      basis_z_eval[ipt + 1*npts] = radial_eval_alpha*x*x*y*z;
      basis_z_eval[ipt + 2*npts] = x*x*(radial_eval + radial_eval_alpha*z*z);
      basis_z_eval[ipt + 3*npts] = radial_eval_alpha*x*y*y*z;
      basis_z_eval[ipt + 4*npts] = x*y*(radial_eval + radial_eval_alpha*z*z);
      basis_z_eval[ipt + 5*npts] = x*z*(2*radial_eval + radial_eval_alpha*z*z);
      basis_z_eval[ipt + 6*npts] = radial_eval_alpha*y*y*y*z;
      basis_z_eval[ipt + 7*npts] = y*y*(radial_eval + radial_eval_alpha*z*z);
      basis_z_eval[ipt + 8*npts] = y*z*(2*radial_eval + radial_eval_alpha*z*z);
      basis_z_eval[ipt + 9*npts] = z*z*(3*radial_eval + radial_eval_alpha*z*z);





#if 0
      // Evaluate the angular part of bfn



      double ang_eval_0;
      double ang_eval_1;
      double ang_eval_2;
      double ang_eval_3;


      ang_eval_0 = radial_eval*x*x*x;
      ang_eval_1 = radial_eval*x*x*y;
      ang_eval_2 = radial_eval*x*x*z;
      ang_eval_3 = radial_eval*x*y*y;
      basis_eval[ipt + 0*npts] = ang_eval_0;
      basis_eval[ipt + 1*npts] = ang_eval_1;
      basis_eval[ipt + 2*npts] = ang_eval_2;
      basis_eval[ipt + 3*npts] = ang_eval_3;

      ang_eval_0 = radial_eval*x*y*z;
      ang_eval_1 = radial_eval*x*z*z;
      ang_eval_2 = radial_eval*y*y*y;
      ang_eval_3 = radial_eval*y*y*z;
      basis_eval[ipt + 4*npts] = ang_eval_0;
      basis_eval[ipt + 5*npts] = ang_eval_1;
      basis_eval[ipt + 6*npts] = ang_eval_2;
      basis_eval[ipt + 7*npts] = ang_eval_3;

      ang_eval_0 = radial_eval*y*z*z;
      ang_eval_1 = radial_eval*z*z*z;
      basis_eval[ipt + 8*npts] = ang_eval_0;
      basis_eval[ipt + 9*npts] = ang_eval_1;


      double dang_eval_x_0, dang_eval_y_0, dang_eval_z_0;
      double dang_eval_x_1, dang_eval_y_1, dang_eval_z_1;
      double dang_eval_x_2, dang_eval_y_2, dang_eval_z_2;
      double dang_eval_x_3, dang_eval_y_3, dang_eval_z_3;

      dang_eval_x_0 = x*x*(3*radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_0 = radial_eval_alpha*x*x*x*y;
      dang_eval_z_0 = radial_eval_alpha*x*x*x*z;
      dang_eval_x_1 = x*y*(2*radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_1 = x*x*(radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_1 = radial_eval_alpha*x*x*y*z;
      dang_eval_x_2 = x*z*(2*radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_2 = radial_eval_alpha*x*x*y*z;
      dang_eval_z_2 = x*x*(radial_eval + radial_eval_alpha*z*z);
      dang_eval_x_3 = y*y*(radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_3 = x*y*(2*radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_3 = radial_eval_alpha*x*y*y*z;
      basis_x_eval[ipt + 0*npts] = dang_eval_x_0;
      basis_y_eval[ipt + 0*npts] = dang_eval_y_0;
      basis_z_eval[ipt + 0*npts] = dang_eval_z_0;
      basis_x_eval[ipt + 1*npts] = dang_eval_x_1;
      basis_y_eval[ipt + 1*npts] = dang_eval_y_1;
      basis_z_eval[ipt + 1*npts] = dang_eval_z_1;
      basis_x_eval[ipt + 2*npts] = dang_eval_x_2;
      basis_y_eval[ipt + 2*npts] = dang_eval_y_2;
      basis_z_eval[ipt + 2*npts] = dang_eval_z_2;
      basis_x_eval[ipt + 3*npts] = dang_eval_x_3;
      basis_y_eval[ipt + 3*npts] = dang_eval_y_3;
      basis_z_eval[ipt + 3*npts] = dang_eval_z_3;

      dang_eval_x_0 = y*z*(radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_0 = x*z*(radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_0 = x*y*(radial_eval + radial_eval_alpha*z*z);
      dang_eval_x_1 = z*z*(radial_eval + radial_eval_alpha*x*x);
      dang_eval_y_1 = radial_eval_alpha*x*y*z*z;
      dang_eval_z_1 = x*z*(2*radial_eval + radial_eval_alpha*z*z);
      dang_eval_x_2 = radial_eval_alpha*x*y*y*y;
      dang_eval_y_2 = y*y*(3*radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_2 = radial_eval_alpha*y*y*y*z;
      dang_eval_x_3 = radial_eval_alpha*x*y*y*z;
      dang_eval_y_3 = y*z*(2*radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_3 = y*y*(radial_eval + radial_eval_alpha*z*z);
      basis_x_eval[ipt + 4*npts] = dang_eval_x_0;
      basis_y_eval[ipt + 4*npts] = dang_eval_y_0;
      basis_z_eval[ipt + 4*npts] = dang_eval_z_0;
      basis_x_eval[ipt + 5*npts] = dang_eval_x_1;
      basis_y_eval[ipt + 5*npts] = dang_eval_y_1;
      basis_z_eval[ipt + 5*npts] = dang_eval_z_1;
      basis_x_eval[ipt + 6*npts] = dang_eval_x_2;
      basis_y_eval[ipt + 6*npts] = dang_eval_y_2;
      basis_z_eval[ipt + 6*npts] = dang_eval_z_2;
      basis_x_eval[ipt + 7*npts] = dang_eval_x_3;
      basis_y_eval[ipt + 7*npts] = dang_eval_y_3;
      basis_z_eval[ipt + 7*npts] = dang_eval_z_3;

      dang_eval_x_0 = radial_eval_alpha*x*y*z*z;
      dang_eval_y_0 = z*z*(radial_eval + radial_eval_alpha*y*y);
      dang_eval_z_0 = y*z*(2*radial_eval + radial_eval_alpha*z*z);
      dang_eval_x_1 = radial_eval_alpha*x*z*z*z;
      dang_eval_y_1 = radial_eval_alpha*y*z*z*z;
      dang_eval_z_1 = z*z*(3*radial_eval + radial_eval_alpha*z*z);
      basis_x_eval[ipt + 8*npts] = dang_eval_x_0;
      basis_y_eval[ipt + 8*npts] = dang_eval_y_0;
      basis_z_eval[ipt + 8*npts] = dang_eval_z_0;
      basis_x_eval[ipt + 9*npts] = dang_eval_x_1;
      basis_y_eval[ipt + 9*npts] = dang_eval_y_1;
      basis_z_eval[ipt + 9*npts] = dang_eval_z_1;

#endif
    } // Loop over points within task
  } // Loop over tasks
        
  } // Loop over shells
} // end kernel

} // namespace GauXC
