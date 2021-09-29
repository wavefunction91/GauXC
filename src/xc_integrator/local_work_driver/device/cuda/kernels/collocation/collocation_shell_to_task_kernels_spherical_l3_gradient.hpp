#pragma once
#include "collocation_device_constants.hpp"
#include "device/xc_device_task.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include "device/common/shell_to_task.hpp"
#include <cassert>

namespace GauXC {


__global__ __launch_bounds__(512,2) void collocation_device_shell_to_task_kernel_spherical_gradient_3(
  uint32_t                         nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks
) {


  __shared__ double alpha[detail::shell_nprim_max], coeff[detail::shell_nprim_max];

  for( auto ish = blockIdx.z; ish < nshell; ish += gridDim.z ) {
  const uint32_t ntasks      = shell_to_task[ish].ntask;
  const auto shell           = shell_to_task[ish].shell_device;
  const auto task_idx        = shell_to_task[ish].task_idx_device;
  const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;


  // Load Shell Data into registers / SM
  const uint32_t nprim = shell->nprim();
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
    for( int ipt = warp_rank /*+ (blockIdx.y*cuda::warp_size)*/; ipt < npts; ipt += (/*gridDim.y * */cuda::warp_size) ) {
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
      double radial_eval_x = 0.;
      double radial_eval_y = 0.;
      double radial_eval_z = 0.;

      #pragma unroll 1
      for( uint32_t i = 0; i < nprim; ++i ) {
        const auto a = alpha[i];
        const auto e = coeff[i] * std::exp( - a * rsq );

        radial_eval += e;

        const auto ae = 2. * a * e;
        radial_eval_x -= ae * x;
        radial_eval_y -= ae * y;
        radial_eval_z -= ae * z;
      }

      
      // Evaluate the angular part of bfn



      double ang_eval_0;
      double ang_eval_1;
      double ang_eval_2;
      double ang_eval_3;


      ang_eval_0 = sqrt_10*radial_eval*y*(3*x*x - y*y)/4;
      ang_eval_1 = sqrt_15*radial_eval*x*y*z;
      ang_eval_2 = sqrt_6*radial_eval*y*(-x*x - y*y + 4*z*z)/4;
      ang_eval_3 = radial_eval*z*(-3*x*x - 3*y*y + 2*z*z)/2;
      basis_eval[ipt + 0*npts] = ang_eval_0;
      basis_eval[ipt + 1*npts] = ang_eval_1;
      basis_eval[ipt + 2*npts] = ang_eval_2;
      basis_eval[ipt + 3*npts] = ang_eval_3;

      ang_eval_0 = sqrt_6*radial_eval*x*(-x*x - y*y + 4*z*z)/4;
      ang_eval_1 = sqrt_15*radial_eval*z*(x*x - y*y)/2;
      ang_eval_2 = sqrt_10*radial_eval*x*(x*x - 3*y*y)/4;
      basis_eval[ipt + 4*npts] = ang_eval_0;
      basis_eval[ipt + 5*npts] = ang_eval_1;
      basis_eval[ipt + 6*npts] = ang_eval_2;


      double dang_eval_x_0, dang_eval_y_0, dang_eval_z_0;
      double dang_eval_x_1, dang_eval_y_1, dang_eval_z_1;
      double dang_eval_x_2, dang_eval_y_2, dang_eval_z_2;
      double dang_eval_x_3, dang_eval_y_3, dang_eval_z_3;

      dang_eval_x_0 = sqrt_10*y*(6*radial_eval*x + radial_eval_x*(3*x*x - y*y))/4;
      dang_eval_y_0 = sqrt_10*(-3*radial_eval*(-x*x + y*y) + radial_eval_y*y*(3*x*x - y*y))/4;
      dang_eval_z_0 = sqrt_10*radial_eval_z*y*(3*x*x - y*y)/4;
      dang_eval_x_1 = sqrt_15*y*z*(radial_eval + radial_eval_x*x);
      dang_eval_y_1 = sqrt_15*x*z*(radial_eval + radial_eval_y*y);
      dang_eval_z_1 = sqrt_15*x*y*(radial_eval + radial_eval_z*z);
      dang_eval_x_2 = -sqrt_6*y*(2*radial_eval*x + radial_eval_x*(x*x + y*y - 4*z*z))/4;
      dang_eval_y_2 = -sqrt_6*(radial_eval*(x*x + 3*y*y - 4*z*z) + radial_eval_y*y*(x*x + y*y - 4*z*z))/4;
      dang_eval_z_2 = sqrt_6*y*(8*radial_eval*z - radial_eval_z*(x*x + y*y - 4*z*z))/4;
      dang_eval_x_3 = -z*(6*radial_eval*x + radial_eval_x*(3*x*x + 3*y*y - 2*z*z))/2;
      dang_eval_y_3 = -z*(6*radial_eval*y + radial_eval_y*(3*x*x + 3*y*y - 2*z*z))/2;
      dang_eval_z_3 = -3*radial_eval*(x*x + y*y - 2*z*z)/2 - radial_eval_z*z*(3*x*x + 3*y*y - 2*z*z)/2;
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

      dang_eval_x_0 = -sqrt_6*(radial_eval*(3*x*x + y*y - 4*z*z) + radial_eval_x*x*(x*x + y*y - 4*z*z))/4;
      dang_eval_y_0 = -sqrt_6*x*(2*radial_eval*y + radial_eval_y*(x*x + y*y - 4*z*z))/4;
      dang_eval_z_0 = sqrt_6*x*(8*radial_eval*z - radial_eval_z*(x*x + y*y - 4*z*z))/4;
      dang_eval_x_1 = sqrt_15*z*(2*radial_eval*x + radial_eval_x*(x*x - y*y))/2;
      dang_eval_y_1 = sqrt_15*z*(-2*radial_eval*y + radial_eval_y*(x*x - y*y))/2;
      dang_eval_z_1 = sqrt_15*(radial_eval + radial_eval_z*z)*(x*x - y*y)/2;
      dang_eval_x_2 = sqrt_10*(3*radial_eval*(x*x - y*y) + radial_eval_x*x*(x*x - 3*y*y))/4;
      dang_eval_y_2 = sqrt_10*x*(-6*radial_eval*y + radial_eval_y*(x*x - 3*y*y))/4;
      dang_eval_z_2 = sqrt_10*radial_eval_z*x*(x*x - 3*y*y)/4;
      basis_x_eval[ipt + 4*npts] = dang_eval_x_0;
      basis_y_eval[ipt + 4*npts] = dang_eval_y_0;
      basis_z_eval[ipt + 4*npts] = dang_eval_z_0;
      basis_x_eval[ipt + 5*npts] = dang_eval_x_1;
      basis_y_eval[ipt + 5*npts] = dang_eval_y_1;
      basis_z_eval[ipt + 5*npts] = dang_eval_z_1;
      basis_x_eval[ipt + 6*npts] = dang_eval_x_2;
      basis_y_eval[ipt + 6*npts] = dang_eval_y_2;
      basis_z_eval[ipt + 6*npts] = dang_eval_z_2;

    } // Loop over points within task
  } // Loop over tasks
  #endif
  } // Loop over shells
} // end kernel

} // namespace GauXC
