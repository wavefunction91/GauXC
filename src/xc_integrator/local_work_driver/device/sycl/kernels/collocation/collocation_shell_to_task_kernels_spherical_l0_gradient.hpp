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
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_vector_types.hpp"
#include "device/common/shell_to_task.hpp"
#include "device/sycl/kernels/sycl_launch.hpp"
#include <sycl/sycl.hpp>
#include <algorithm>
#include <cassert>

namespace GauXC {

// Launch geometry inherited from the CUDA backend: 512 work-items per
// work-group, 16 sub-groups worth of shell scratch in local memory
static constexpr uint32_t collocation_device_shell_to_task_kernel_spherical_gradient_0_max_wg = 512;
static constexpr uint32_t collocation_device_shell_to_task_kernel_spherical_gradient_0_nwarp  = 16;


void collocation_device_shell_to_task_kernel_spherical_gradient_0(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {

  // Recover the work-item and shell scratch the way the CUDA kernel
  // reads threadIdx and declares __shared__.
  auto it = GauXC::sycl::this_item();
  auto& alpha = GauXC::sycl::local_mem<double[16][detail::shell_nprim_max + 1]>();
  auto& coeff = GauXC::sycl::local_mem<double[16][detail::shell_nprim_max + 1]>();

  // Shell scratch lives in work-group local memory, one row per sub-group
  const uint32_t local_warp_id = it.get_local_id(2) / sycl::warp_size;
  double* my_alpha = &alpha[local_warp_id][0];
  double* my_coeff = &coeff[local_warp_id][0];

  for( auto ish = it.get_group(0); ish < nshell; ish += it.get_group_range(0) ) {
  const uint32_t ntasks      = shell_to_task[ish].ntask;
  const auto shell           = shell_to_task[ish].shell_device;
  const auto task_idx        = shell_to_task[ish].task_idx_device;
  const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;


  // Load Shell Data into registers / SM
  const uint32_t nprim = shell->nprim();
  const double3 O  = *reinterpret_cast<const double3*>(shell->O_data());

  const int global_warp_id = (it.get_global_id(2)) / sycl::warp_size;
  const int nwarp_global   = std::max<int>(it.get_global_range(2) / sycl::warp_size,1);

  // Read in coeffs/exps into SM on first warp
  {
    auto* coeff_gm = shell->coeff_data();
    auto* alpha_gm = shell->alpha_data();
    static_assert( detail::shell_nprim_max == sycl::warp_size );
    const int warp_rank = it.get_local_id(2) % sycl::warp_size;
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
        for( int ipt = it.get_local_id(2) % sycl::warp_size; ipt < npts; ipt += sycl::warp_size ) {
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

            for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = my_alpha[i];
        const auto e = my_coeff[i] * ::sycl::exp( - a * rsq );

        radial_eval += e;
        radial_eval_alpha += a * e;
      }

      radial_eval_alpha *= -2;

      // Common Subexpressions


      // Evaluate basis function
      basis_eval[ipt + 0*npts] = radial_eval;


    
      // Evaluate first derivative of bfn wrt x
      basis_x_eval[ipt + 0*npts] = radial_eval_alpha*x;

      // Evaluate first derivative of bfn wrt y
      basis_y_eval[ipt + 0*npts] = radial_eval_alpha*y;

      // Evaluate first derivative of bfn wrt z
      basis_z_eval[ipt + 0*npts] = radial_eval_alpha*z;







#if 0
      // Evaluate the angular part of bfn



      double ang_eval_0;


      ang_eval_0 = radial_eval;
      basis_eval[ipt + 0*npts] = ang_eval_0;


      double dang_eval_x_0, dang_eval_y_0, dang_eval_z_0;

      dang_eval_x_0 = radial_eval_alpha*x;
      dang_eval_y_0 = radial_eval_alpha*y;
      dang_eval_z_0 = radial_eval_alpha*z;
      basis_x_eval[ipt + 0*npts] = dang_eval_x_0;
      basis_y_eval[ipt + 0*npts] = dang_eval_y_0;
      basis_z_eval[ipt + 0*npts] = dang_eval_z_0;

#endif
    } // Loop over points within task
  } // Loop over tasks
        
  } // Loop over shells
} // end kernel

} // namespace GauXC
