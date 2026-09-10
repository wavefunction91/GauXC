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

static constexpr uint32_t collocation_device_shell_to_task_kernel_cartesian_2_max_wg = 512;
static constexpr uint32_t collocation_device_shell_to_task_kernel_cartesian_2_nwarp  = 16;


void collocation_device_shell_to_task_kernel_cartesian_2(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {

  auto it = GauXC::sycl::this_item();
  auto& alpha = GauXC::sycl::local_mem<double[16][detail::shell_nprim_max + 1]>();
  auto& coeff = GauXC::sycl::local_mem<double[16][detail::shell_nprim_max + 1]>();

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

            for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = my_alpha[i];
        const auto e = my_coeff[i] * ::sycl::exp( - a * rsq );

        radial_eval += e;
      }


      // Common Subexpressions
      const auto x0 = radial_eval*x; 


      // Evaluate basis function
      basis_eval[ipt + 0*npts] = radial_eval*(x*x);
      basis_eval[ipt + 1*npts] = x0*y;
      basis_eval[ipt + 2*npts] = x0*z;
      basis_eval[ipt + 3*npts] = radial_eval*(y*y);
      basis_eval[ipt + 4*npts] = radial_eval*y*z;
      basis_eval[ipt + 5*npts] = radial_eval*(z*z);


    







#if 0
      // Evaluate the angular part of bfn



      double ang_eval_0;
      double ang_eval_1;
      double ang_eval_2;
      double ang_eval_3;


      ang_eval_0 = radial_eval*(x*x);
      ang_eval_1 = x0*y;
      ang_eval_2 = x0*z;
      ang_eval_3 = radial_eval*(y*y);
      basis_eval[ipt + 0*npts] = ang_eval_0;
      basis_eval[ipt + 1*npts] = ang_eval_1;
      basis_eval[ipt + 2*npts] = ang_eval_2;
      basis_eval[ipt + 3*npts] = ang_eval_3;

      ang_eval_0 = radial_eval*y*z;
      ang_eval_1 = radial_eval*(z*z);
      basis_eval[ipt + 4*npts] = ang_eval_0;
      basis_eval[ipt + 5*npts] = ang_eval_1;


#endif
    } // Loop over points within task
  } // Loop over tasks
        
  } // Loop over shells
} // end kernel

} // namespace GauXC
