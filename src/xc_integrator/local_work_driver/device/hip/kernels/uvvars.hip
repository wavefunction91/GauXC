/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/uvvars.hpp"
#include "hip_extensions.hpp"
#include "device_specific/hip_device_constants.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/hip_util.hpp"

namespace GauXC {


__global__ void eval_uvars_lda_kernel( size_t           ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  auto* den_eval_device   = task.den;

  const auto* basis_eval_device = task.bf;

  const auto* den_basis_prod_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  double den_reg = 0.;

  if( tid_x < nbf and tid_y < npts ) {

    const double* bf_col   = basis_eval_device     + tid_x*npts;
    const double* db_col   = den_basis_prod_device + tid_x*npts;

    den_reg = bf_col[ tid_y ]   * db_col[ tid_y ];

  }

  // Warp blocks are stored col major
  den_reg = 2 * hip::warp_reduce_sum<hip::warp_size>( den_reg );


  if( threadIdx.x == 0 and tid_y < npts ) {
    atomicAdd( den_eval_device   + tid_y, den_reg );
  }


}


void eval_uvars_lda( size_t ntasks, int32_t npts_max,
                   integrator_ks_scheme scheme,
  XCDeviceTask* device_tasks,
 device_queue queue ) {
  hipStream_t stream = queue.queue_as<util::hip_stream>();

  dim3 threads(hip::warp_size, hip::max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( npts_max , threads.x ),
               util::div_ceil( npts_max , threads.y ),
               ntasks );

  hipLaunchKernelGGL(eval_uvars_lda_kernel, dim3(blocks), dim3(threads), 0, 
    stream, ntasks, device_tasks );

}













__global__ void eval_uvars_gga_kernel( size_t           ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  auto* den_eval_device   = task.den;
  auto* den_x_eval_device = task.den_x;
  auto* den_y_eval_device = task.den_y;
  auto* den_z_eval_device = task.den_z;

  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  const auto* den_basis_prod_device = task.zmat;

  // We always launch enough blocks to cover npts, so blocks aren't doing multiple results
  double den_reg = 0.;
  double dx_reg = 0.;
  double dy_reg = 0.;
  double dz_reg = 0.;

  // Have each thread accumulate its own reduction result into a register.
  // There's no real _need_ for LDS because the reductions are small and
  // therefore can be done without sharing.
  for( int ibf = 0; ibf < nbf; ibf++ ) {

    for( int  ipt = blockIdx.x * blockDim.x + threadIdx.x; ipt < npts; ipt += blockDim.x * gridDim.x ) {

      const double* bf_col   = basis_eval_device     + ibf*npts;
      const double* bf_x_col = dbasis_x_eval_device  + ibf*npts;
      const double* bf_y_col = dbasis_y_eval_device  + ibf*npts;
      const double* bf_z_col = dbasis_z_eval_device  + ibf*npts;
      const double* db_col   = den_basis_prod_device + ibf*npts;

      den_reg += 2 * bf_col[ ipt ]   * db_col[ ipt ];
      dx_reg += 4 * bf_x_col[ ipt ] * db_col[ ipt ];
      dy_reg += 4 * bf_y_col[ ipt ] * db_col[ ipt ];
      dz_reg += 4 * bf_z_col[ ipt ] * db_col[ ipt ];
    }
  }


  for( int  ipt = blockIdx.x * blockDim.x + threadIdx.x; ipt < npts; ipt += blockDim.x * gridDim.x ) {
    den_eval_device   [ipt] = den_reg;
    den_x_eval_device [ipt] = dx_reg ;
    den_y_eval_device [ipt] = dy_reg ;
    den_z_eval_device [ipt] = dz_reg ;
  }

}


__global__ void eval_vvars_gga_kernel( 
  size_t        npts,
  const double* den_x_eval_device,
  const double* den_y_eval_device,
  const double* den_z_eval_device,
        double* gamma_eval_device
) {

  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < npts ) {

    const double dx = den_x_eval_device[ tid ];
    const double dy = den_y_eval_device[ tid ];
    const double dz = den_z_eval_device[ tid ];

    gamma_eval_device[tid] = dx*dx + dy*dy + dz*dz;

  }

}




void eval_uvvars_gga( size_t ntasks, size_t npts_total, int32_t nbf_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>();

  // U Variables
  {
  dim3 threads(hip::max_threads_per_thread_block, 1, 1);
  dim3 blocks( util::div_ceil( npts_max , threads.x ),
               1,
               ntasks );

  hipLaunchKernelGGL(eval_uvars_gga_kernel, dim3(blocks), dim3(threads), 0, 
    stream, ntasks, device_tasks );
  }

  // V Variables
  dim3 threads( hip::max_threads_per_thread_block );
  dim3 blocks( util::div_ceil( npts_total, threads.x ) );
  hipLaunchKernelGGL(eval_vvars_gga_kernel, blocks, threads, 0, stream,
    npts_total, denx, deny, denz, gamma);

}

#define GGA_KERNEL_SM_BLOCK_Y 16

template <bool need_lapl>
__global__ void eval_uvars_mgga_kernel( size_t           ntasks,
                                       XCDeviceTask* tasks_device ) {

  constexpr auto warp_size = hip::warp_size;
  //constexpr auto max_warps_per_thread_block = hip::max_warps_per_thread_block;

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  auto* tau_eval_device   = task.tau;
  decltype(tau_eval_device) lapl_eval_device = nullptr;
  if constexpr (need_lapl) {
    lapl_eval_device = task.denlapl;
  }

  //const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;
  decltype(dbasis_x_eval_device) basis_lapl_eval_device = nullptr;
  if constexpr (need_lapl) {
    basis_lapl_eval_device = task.d2bflapl;
  }

  //const auto* den_basis_prod_device    = task.zmat;
  const auto* den_basis_dx_prod_device = task.xmat_x;
  const auto* den_basis_dy_prod_device = task.xmat_y;
  const auto* den_basis_dz_prod_device = task.xmat_z;
  decltype(den_basis_dx_prod_device) den_basis_prod_device = nullptr;
  if constexpr (need_lapl) {
    den_basis_prod_device = task.zmat;
  }

  __shared__ double den_shared[3+!!need_lapl][warp_size][GGA_KERNEL_SM_BLOCK_Y+1];

  for ( int bid_x = blockIdx.x * blockDim.x; 
        bid_x < nbf;
        bid_x += blockDim.x * gridDim.x ) {
    
    for ( int bid_y = blockIdx.y * GGA_KERNEL_SM_BLOCK_Y; 
          bid_y < npts;
          bid_y += GGA_KERNEL_SM_BLOCK_Y * gridDim.y ) {
        
      for (int sm_y = threadIdx.y; sm_y < GGA_KERNEL_SM_BLOCK_Y; sm_y += blockDim.y) {
        den_shared[0][threadIdx.x][sm_y] = 0.;
        den_shared[1][threadIdx.x][sm_y] = 0.;
        den_shared[2][threadIdx.x][sm_y] = 0.;
        if constexpr (need_lapl)
          den_shared[3][threadIdx.x][sm_y] = 0.;

        if (bid_y + threadIdx.x < npts and bid_x + sm_y < nbf) { 
          const double* db_x_col = den_basis_dx_prod_device + (bid_x + sm_y)*npts;
          const double* db_y_col = den_basis_dy_prod_device + (bid_x + sm_y)*npts;
          const double* db_z_col = den_basis_dz_prod_device + (bid_x + sm_y)*npts;

          const double* bf_x_col = dbasis_x_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_y_col = dbasis_y_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_z_col = dbasis_z_eval_device  + (bid_x + sm_y)*npts;


          den_shared[0][threadIdx.x][sm_y] = bf_x_col[ bid_y + threadIdx.x ] * db_x_col[ bid_y + threadIdx.x ];
          den_shared[1][threadIdx.x][sm_y] = bf_y_col[ bid_y + threadIdx.x ] * db_y_col[ bid_y + threadIdx.x ];
          den_shared[2][threadIdx.x][sm_y] = bf_z_col[ bid_y + threadIdx.x ] * db_z_col[ bid_y + threadIdx.x ];


          if constexpr (need_lapl) {
            const double* db_col   = den_basis_prod_device  + (bid_x + sm_y)*npts;
            const double* bf_l_col = basis_lapl_eval_device + (bid_x + sm_y)*npts;
            den_shared[3][threadIdx.x][sm_y] = bf_l_col[ bid_y + threadIdx.x ] * db_col[ bid_y + threadIdx.x ];
          }
        }
      }
      __syncthreads();


      for (int sm_y = threadIdx.y; sm_y < GGA_KERNEL_SM_BLOCK_Y; sm_y += blockDim.y) {
        const int tid_y = bid_y + sm_y;

        double tx_reg  = den_shared[0][sm_y][threadIdx.x];
        double ty_reg  = den_shared[1][sm_y][threadIdx.x];
        double tz_reg  = den_shared[2][sm_y][threadIdx.x];
        // Warp blocks are stored col major
        double tau_reg = 0.0;
        tau_reg  = 0.5 * hip::warp_reduce_sum<warp_size>( tx_reg );
        tau_reg += 0.5 * hip::warp_reduce_sum<warp_size>( ty_reg );
        tau_reg += 0.5 * hip::warp_reduce_sum<warp_size>( tz_reg );

        double lapl_reg = 0.0;
        if constexpr (need_lapl) {
          lapl_reg = den_shared[3][sm_y][threadIdx.x];
          lapl_reg = hip::warp_reduce_sum<warp_size>(lapl_reg);
          lapl_reg = 2. * lapl_reg + 4. * tau_reg;
        }

        if( threadIdx.x == 0 and tid_y < npts ) {
          atomicAdd( tau_eval_device   + tid_y, tau_reg );
          if constexpr (need_lapl) {
            atomicAdd( lapl_eval_device   + tid_y, lapl_reg );
          }
        }
      }
      __syncthreads();
    }
  }
}


void eval_uvvars_mgga( size_t ntasks, size_t npts_total, int32_t nbf_max, 
  int32_t npts_max, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, bool do_lapl,
  device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>();

  // U Variables
  {
  dim3 threads( hip::warp_size, hip::max_warps_per_thread_block / 2, 1 );
  dim3 blocks( std::min(uint64_t(4), util::div_ceil( nbf_max, 4 )),
               std::min(uint64_t(GGA_KERNEL_SM_BLOCK_Y), util::div_ceil( npts_max, GGA_KERNEL_SM_BLOCK_Y )),
               ntasks );
  eval_uvars_gga_kernel <<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
  if(do_lapl)
    eval_uvars_mgga_kernel<true><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
  else
    eval_uvars_mgga_kernel<false><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
  }

  // V variables (GAMMA)
  dim3 threads( hip::max_threads_per_thread_block );
  dim3 blocks( util::div_ceil( npts_total, threads.x ) );
  eval_vvars_gga_kernel<<< blocks, threads, 0, stream >>>(
    npts_total, denx, deny, denz, gamma
  );
}


}
