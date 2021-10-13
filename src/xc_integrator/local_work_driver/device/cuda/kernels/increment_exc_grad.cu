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

      const auto* __restrict__  zmat = task->zmat + shoff;
      const auto* __restrict__  vrho = task->vrho;

      #pragma unroll 1
      for( uint32_t ipt = threadIdx.x % cuda::warp_size; 
           ipt < npts; 
           ipt += cuda::warp_size ) {

        const double vrho_i = vrho[ipt];
        for( uint32_t ibf = 0; ibf < shsz; ++ibf ) {
          const double z_mu_i    = vrho_i * 2. * zmat[ipt + ibf*npts];
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


}
