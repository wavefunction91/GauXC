/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "device_specific/cuda_device_constants.hpp"
#include "device_specific/cuda_util.hpp"
#include "device/xc_device_data.hpp"

#define MGGA_KERNEL_SM_BLOCK 32

namespace GauXC {



template <density_id den_select, bool need_lapl>
__global__ void eval_vvar_mgga_kern( size_t           ntasks,
                                     XCDeviceTask* tasks_device ) {

  constexpr auto warp_size = cuda::warp_size;
  //constexpr auto max_warps_per_thread_block = cuda::max_warps_per_thread_block;

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  double* tau_eval_device  = nullptr;
  double* lapl_eval_device = nullptr;

  if constexpr (den_select == DEN_S) {
    tau_eval_device = task.tau_s;
    if constexpr (need_lapl) {
      lapl_eval_device = task.lapl_s;
    }
  }

  if constexpr (den_select == DEN_Z) {
    tau_eval_device = task.tau_z;
    if constexpr (need_lapl) {
      lapl_eval_device = task.lapl_z;
    }
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

  __shared__ double den_shared[3+!!need_lapl][warp_size][MGGA_KERNEL_SM_BLOCK+1];

  for ( int bid_x = blockIdx.x * blockDim.x; 
        bid_x < nbf;
        bid_x += blockDim.x * gridDim.x ) {
    
    for ( int bid_y = blockIdx.y * MGGA_KERNEL_SM_BLOCK; 
          bid_y < npts;
          bid_y += MGGA_KERNEL_SM_BLOCK * gridDim.y ) {
        
      for (int sm_y = threadIdx.y; sm_y < MGGA_KERNEL_SM_BLOCK; sm_y += blockDim.y) {
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


      for (int sm_y = threadIdx.y; sm_y < MGGA_KERNEL_SM_BLOCK; sm_y += blockDim.y) {
        const int tid_y = bid_y + sm_y;

        register double tx_reg  = den_shared[0][sm_y][threadIdx.x];
        register double ty_reg  = den_shared[1][sm_y][threadIdx.x];
        register double tz_reg  = den_shared[2][sm_y][threadIdx.x];
        // Warp blocks are stored col major
        register double tau_reg = 0.0;
        tau_reg  = 0.5 * cuda::warp_reduce_sum<warp_size>( tx_reg );
        tau_reg += 0.5 * cuda::warp_reduce_sum<warp_size>( ty_reg );
        tau_reg += 0.5 * cuda::warp_reduce_sum<warp_size>( tz_reg );

        register double lapl_reg = 0.0;
        if constexpr (need_lapl) {
          lapl_reg = den_shared[3][sm_y][threadIdx.x];
          lapl_reg = cuda::warp_reduce_sum<warp_size>(lapl_reg);
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




template <bool need_lapl>
__global__ void eval_uvars_mgga_uks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts  = task.npts;

  auto* tau_pos_eval_device = task.tau_s;
  auto* tau_neg_eval_device = task.tau_z;

  double* lapl_pos_eval_device = nullptr;
  double* lapl_neg_eval_device = nullptr;
  if constexpr (need_lapl) {
    lapl_pos_eval_device = task.lapl_s;
    lapl_neg_eval_device = task.lapl_z;
  }

  const int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if( tid < npts ) {
    const double ts = tau_pos_eval_device[ tid ];
    const double tz = tau_neg_eval_device[ tid ];
    tau_pos_eval_device[ tid ] = 0.5*(ts + tz);
    tau_neg_eval_device[ tid ] = 0.5*(ts - tz);

    if constexpr (need_lapl) {
      const double ls = lapl_pos_eval_device[ tid ];
      const double lz = lapl_neg_eval_device[ tid ];
      lapl_pos_eval_device[ tid ] = 0.5*(ls + lz);
      lapl_neg_eval_device[ tid ] = 0.5*(ls - lz);
    }
  }

}

}
