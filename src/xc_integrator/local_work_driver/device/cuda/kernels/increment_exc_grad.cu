/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/increment_exc_grad.hpp"
#include "cuda_extensions.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"

namespace GauXC {

__global__ __launch_bounds__(1024,1) void increment_exc_grad_lda_kernel(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks,
  double*            __restrict__ EXC_GRAD
) {

  for( uint32_t ish = blockIdx.z; ish < nshell; ish += gridDim.z ) {
    const uint32_t ntasks      = shell_to_task[ish].ntask;
    const auto shell           = shell_to_task[ish].shell_device;
    const auto task_idx        = shell_to_task[ish].task_idx_device;
    const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;
    const uint32_t shsz   = shell->size();

    const int global_warp_id = 
      (threadIdx.x + blockIdx.x*blockDim.x) / cuda::warp_size;
    const int nwarp_global   = max((blockDim.x*gridDim.x) / cuda::warp_size,1);

    double g_acc_x(0), g_acc_y(0), g_acc_z(0);
    for( uint32_t itask = global_warp_id; itask < ntasks; itask += nwarp_global ) {
      
      const auto*    task   = device_tasks + task_idx[itask];
      const uint32_t npts   = task->npts;
      const size_t   shoff  = task_shell_offs[itask] * npts;

      const auto* __restrict__ basis_x_eval = task->dbfx + shoff;
      const auto* __restrict__ basis_y_eval = task->dbfy + shoff;
      const auto* __restrict__ basis_z_eval = task->dbfz + shoff;

      const auto* __restrict__  xmat = task->zmat + shoff;
      const auto* __restrict__  vrho = task->vrho;

      #pragma unroll 1
      for( uint32_t ipt = threadIdx.x % cuda::warp_size; 
           ipt < npts; 
           ipt += cuda::warp_size ) {

        const double vrho_i = vrho[ipt];
        for( uint32_t ibf = 0; ibf < shsz; ++ibf ) {
          const double z_mu_i    = vrho_i * xmat[ipt + ibf*npts];
          const double dbfx_mu_i = basis_x_eval[ipt + ibf*npts];
          const double dbfy_mu_i = basis_y_eval[ipt + ibf*npts];
          const double dbfz_mu_i = basis_z_eval[ipt + ibf*npts];

          g_acc_x += z_mu_i * dbfx_mu_i;
          g_acc_y += z_mu_i * dbfy_mu_i;
          g_acc_z += z_mu_i * dbfz_mu_i;
        } // Loop over bfns within a shell

      } // Loop over points

    } // Loop over tasks assigned to shell

    constexpr auto warp_size = cuda::warp_size;
    g_acc_x = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_x );
    g_acc_y = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_y );
    g_acc_z = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_z );

    if( (threadIdx.x % cuda::warp_size) == 0 ) {
      const int iCen = shell_to_task[ish].center_idx;
      
      atomicAdd( EXC_GRAD + 3*iCen + 0, g_acc_x );
      atomicAdd( EXC_GRAD + 3*iCen + 1, g_acc_y );
      atomicAdd( EXC_GRAD + 3*iCen + 2, g_acc_z );
    }

  } // Loop over shells

}

void increment_exc_grad_lda( size_t nshell, ShellToTaskDevice* shell_to_task,
  XCDeviceTask* device_tasks, double* EXC_GRAD, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  #if 0
  int nthreads_per_block = 1024;
  int nwarp_per_block    = nthreads_per_block / cuda::warp_size;
  int nblocks            = util::div_ceil( nshell, nwarp_per_block );

  dim3 threads( nthreads_per_block );
  dim3 blocks( nblocks );
  #else
  dim3 threads(1024), blocks(1,1,nshell);
  #endif

  increment_exc_grad_lda_kernel<<<blocks, threads, 0 , stream>>>(
    nshell, shell_to_task, device_tasks, EXC_GRAD 
  );
}



















__global__ __launch_bounds__(512,1) void increment_exc_grad_gga_kernel(
  uint32_t                        nshell,
  ShellToTaskDevice* __restrict__ shell_to_task,
  XCDeviceTask*      __restrict__ device_tasks,
  double*            __restrict__ EXC_GRAD
) {

  for( uint32_t ish = blockIdx.z; ish < nshell; ish += gridDim.z ) {
    const uint32_t ntasks      = shell_to_task[ish].ntask;
    const auto shell           = shell_to_task[ish].shell_device;
    const auto task_idx        = shell_to_task[ish].task_idx_device;
    const auto task_shell_offs = shell_to_task[ish].task_shell_offs_device;
    const uint32_t shsz   = shell->size();

    const int global_warp_id = 
      (threadIdx.x + blockIdx.x*blockDim.x) / cuda::warp_size;
    const int nwarp_global   = max((blockDim.x*gridDim.x) / cuda::warp_size,1);

    double g_acc_x(0), g_acc_y(0), g_acc_z(0);
    for( uint32_t itask = global_warp_id; itask < ntasks; itask += nwarp_global ) {
      
      const auto*    task   = device_tasks + task_idx[itask];
      const uint32_t npts   = task->npts;
      const size_t   shoff  = task_shell_offs[itask] * npts;

      const auto* __restrict__ basis_x_eval = task->dbfx + shoff;
      const auto* __restrict__ basis_y_eval = task->dbfy + shoff;
      const auto* __restrict__ basis_z_eval = task->dbfz + shoff;

      const auto* __restrict__ basis_xx_eval = task->d2bfxx + shoff;
      const auto* __restrict__ basis_xy_eval = task->d2bfxy + shoff;
      const auto* __restrict__ basis_xz_eval = task->d2bfxz + shoff;
      const auto* __restrict__ basis_yy_eval = task->d2bfyy + shoff;
      const auto* __restrict__ basis_yz_eval = task->d2bfyz + shoff;
      const auto* __restrict__ basis_zz_eval = task->d2bfzz + shoff;

      const auto* __restrict__  xmat = task->zmat + shoff;
      const auto* __restrict__  xmat_x = task->xmat_x + shoff;
      const auto* __restrict__  xmat_y = task->xmat_y + shoff;
      const auto* __restrict__  xmat_z = task->xmat_z + shoff;

      const auto* __restrict__  vrho = task->vrho;
      const auto* __restrict__  vgamma = task->vgamma;

      const auto* __restrict__ den_x = task->ddenx;
      const auto* __restrict__ den_y = task->ddeny;
      const auto* __restrict__ den_z = task->ddenz;

      #pragma unroll 1
      for( uint32_t ipt = threadIdx.x % cuda::warp_size; 
           ipt < npts; 
           ipt += cuda::warp_size ) {

        const double vrho_i   = vrho[ipt];
        const double vgamma_i = vgamma[ipt];

        const double denx_i = den_x[ipt];
        const double deny_i = den_y[ipt];
        const double denz_i = den_z[ipt];
        for( uint32_t ibf = 0; ibf < shsz; ++ibf ) {
          const double z_mu_i    = xmat[ipt + ibf*npts];
          const double dbfx_mu_i = basis_x_eval[ipt + ibf*npts];
          const double dbfy_mu_i = basis_y_eval[ipt + ibf*npts];
          const double dbfz_mu_i = basis_z_eval[ipt + ibf*npts];

          g_acc_x += vrho_i * z_mu_i * dbfx_mu_i;
          g_acc_y += vrho_i * z_mu_i * dbfy_mu_i;
          g_acc_z += vrho_i * z_mu_i * dbfz_mu_i;

          const double zx = xmat_x[ipt + ibf*npts];
          const double zy = xmat_y[ipt + ibf*npts];
          const double zz = xmat_z[ipt + ibf*npts];

	        const double d11_xmat_term = denx_i * zx + deny_i * zy + denz_i * zz;

          const double d2bfxx = basis_xx_eval[ipt + ibf*npts];
          const double d2bfxy = basis_xy_eval[ipt + ibf*npts];
          const double d2bfxz = basis_xz_eval[ipt + ibf*npts];
          const double d2bfyy = basis_yy_eval[ipt + ibf*npts];
          const double d2bfyz = basis_yz_eval[ipt + ibf*npts];
          const double d2bfzz = basis_zz_eval[ipt + ibf*npts];

	        const double d2_term_x = d2bfxx*denx_i + d2bfxy*deny_i + d2bfxz*denz_i;
	        const double d2_term_y = d2bfxy*denx_i + d2bfyy*deny_i + d2bfyz*denz_i;
	        const double d2_term_z = d2bfxz*denx_i + d2bfyz*deny_i + d2bfzz*denz_i;

	        g_acc_x += 2 * vgamma_i * ( z_mu_i * d2_term_x + dbfx_mu_i * d11_xmat_term );
	        g_acc_y += 2 * vgamma_i * ( z_mu_i * d2_term_y + dbfy_mu_i * d11_xmat_term );
	        g_acc_z += 2 * vgamma_i * ( z_mu_i * d2_term_z + dbfz_mu_i * d11_xmat_term );

        } // Loop over bfns within a shell

      } // Loop over points

    } // Loop over tasks assigned to shell

    constexpr auto warp_size = cuda::warp_size;
    g_acc_x = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_x );
    g_acc_y = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_y );
    g_acc_z = -2. * cuda::warp_reduce_sum<warp_size>( g_acc_z );

    if( (threadIdx.x % cuda::warp_size) == 0 ) {
      const int iCen = shell_to_task[ish].center_idx;
      
      atomicAdd( EXC_GRAD + 3*iCen + 0, g_acc_x );
      atomicAdd( EXC_GRAD + 3*iCen + 1, g_acc_y );
      atomicAdd( EXC_GRAD + 3*iCen + 2, g_acc_z );
    }

  } // Loop over shells

}

void increment_exc_grad_gga( size_t nshell, ShellToTaskDevice* shell_to_task,
  XCDeviceTask* device_tasks, double* EXC_GRAD, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads(512), blocks(1,1,nshell);

  increment_exc_grad_gga_kernel<<<blocks, threads, 0 , stream>>>(
    nshell, shell_to_task, device_tasks, EXC_GRAD 
  );
}

}
