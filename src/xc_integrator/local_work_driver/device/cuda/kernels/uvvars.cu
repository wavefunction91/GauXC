/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/uvvars.hpp"
#include "cuda_extensions.hpp"
#include "device_specific/cuda_device_constants.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"
#include "device/xc_device_data.hpp"

namespace GauXC {

#define GGA_KERNEL_SM_BLOCK_Y 32

__global__ void eval_uvars_lda_rks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {
  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;
  
  const auto& task = tasks_device[ batch_idx ];
  const auto npts  = task.npts;
  
  auto*         den_eval_device     = task.den;
  auto*         den_s_eval_device   = task.den_s;

  const int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if( tid < npts ) {
    const double ps = den_s_eval_device[ tid ];
    den_eval_device  [ tid ] = ps;

  }

}

__global__ void eval_uvars_lda_uks_kernel( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;

  auto* den_pos_eval_device   = task.den_s;
  auto* den_neg_eval_device   = task.den_z;


  const int tid_y = blockIdx.x * blockDim.x + threadIdx.x;


  if( tid_y < npts ) {
    const auto ps = den_pos_eval_device[ tid_y ];
    const auto pz = den_neg_eval_device[ tid_y ];
    den_pos_eval_device[ tid_y ] = 0.5*(ps + pz);
    den_neg_eval_device[ tid_y ] = 0.5*(ps - pz);

  }
}

__global__ void eval_uvars_lda_gks_kernel( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;

  auto* den_z_eval_device   = task.den_s;
  auto* den_s_eval_device   = task.den_z;
  auto* den_y_eval_device   = task.den_y;
  auto* den_x_eval_device   = task.den_x;
  auto* K_z_eval_device     = task.K_z;
  auto* K_y_eval_device     = task.K_y;
  auto* K_x_eval_device     = task.K_x;
  const double dtolsq = 1e-24;  // TODO: make variable

  const int tid = blockIdx.x * blockDim.x + threadIdx.x;


  if( tid < npts ) {
    const auto ps = den_s_eval_device[ tid ];
    const auto pz = den_z_eval_device[ tid ];
    const auto py = den_y_eval_device[ tid ];
    const auto px = den_x_eval_device[ tid ];
    const auto mtemp = pz*pz + px*px + py*py;
    double mnorm = 0.;
  
    if (mtemp > dtolsq) {
      mnorm = sqrt(mtemp);
      K_z_eval_device[ tid ] = pz / mnorm;
      K_y_eval_device[ tid ] = py / mnorm;
      K_x_eval_device[ tid ] = px / mnorm;
    }
    else {
      mnorm = (1. / 3.) * (px + py + pz);
      K_z_eval_device[ tid ] = 1. / 3.;
      K_y_eval_device[ tid ] = 1. / 3.;
      K_x_eval_device[ tid ] = 1. / 3.;
    }


    den_s_eval_device[ tid ] = 0.5*(ps + mnorm);
    den_z_eval_device[ tid ] = 0.5*(ps - mnorm);

  }
}


__global__ void eval_uvars_gga_rks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {
  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;
  
  const auto& task = tasks_device[ batch_idx ];
  const auto npts  = task.npts;
  
  const auto*   dden_sx_eval_device = task.dden_sx;
  const auto*   dden_sy_eval_device = task.dden_sy;
  const auto*   dden_sz_eval_device = task.dden_sz;
  auto*         gamma_eval_device   = task.gamma;

  const int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if( tid < npts ) {
    const double dx = dden_sx_eval_device[ tid ];
    const double dy = dden_sy_eval_device[ tid ];
    const double dz = dden_sz_eval_device[ tid ];

    gamma_eval_device[ tid ] = dx*dx + dy*dy + dz*dz;

  }

}

__global__ void eval_uvars_gga_uks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;

  auto*           den_pos_eval_device   = task.den_s;
  const auto*     den_pos_x_eval_device = task.dden_sx;
  const auto*     den_pos_y_eval_device = task.dden_sy;
  const auto*     den_pos_z_eval_device = task.dden_sz;

  auto*           den_neg_eval_device   = task.den_z;
  const auto*     den_neg_x_eval_device = task.dden_zx;
  const auto*     den_neg_y_eval_device = task.dden_zy;
  const auto*     den_neg_z_eval_device = task.dden_zz;

  auto*     gamma_pp_eval_device  = task.gamma_pp;
  auto*     gamma_pm_eval_device  = task.gamma_pm;
  auto*     gamma_mm_eval_device  = task.gamma_mm;

  const int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if( tid < npts ) {
    const double ps     = den_pos_eval_device[ tid ];
    const double pz     = den_neg_eval_device[ tid ];
    const double dndx   = den_pos_x_eval_device[ tid ];
    const double dndy   = den_pos_y_eval_device[ tid ];
    const double dndz   = den_pos_z_eval_device[ tid ];
    const double dMzdx  = den_neg_x_eval_device[ tid ];
    const double dMzdy  = den_neg_y_eval_device[ tid ];
    const double dMzdz  = den_neg_z_eval_device[ tid ];

    // (del n).(del n)
    const auto dn_sq  = dndx*dndx + dndy*dndy + dndz*dndz;
    // (del Mz).(del Mz)
    const auto dMz_sq = dMzdx*dMzdx + dMzdy*dMzdy + dMzdz*dMzdz;
    // (del n).(del Mz)
    const auto dn_dMz = dndx*dMzdx + dndy*dMzdy + dndz*dMzdz;

    gamma_pp_eval_device[ tid ] = 0.25*(dn_sq + dMz_sq) + 0.5*dn_dMz;
    gamma_pm_eval_device[ tid ] = 0.25*(dn_sq - dMz_sq);
    gamma_mm_eval_device[ tid ] = 0.25*(dn_sq + dMz_sq) - 0.5*dn_dMz;

    den_pos_eval_device[ tid ] = 0.5*(ps + pz);
    den_neg_eval_device[ tid ] = 0.5*(ps - pz);
  }

}

__global__ void eval_uvars_gga_gks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;

        auto*     den_s_eval_device   = task.den_s;
  const auto*     dden_sx_eval_device = task.dden_sx;
  const auto*     dden_sy_eval_device = task.dden_sy;
  const auto*     dden_sz_eval_device = task.dden_sz;

        auto*     den_z_eval_device   = task.den_z;
  const auto*     dden_zx_eval_device = task.dden_zx;
  const auto*     dden_zy_eval_device = task.dden_zy;
  const auto*     dden_zz_eval_device = task.dden_zz;

  const auto*     den_y_eval_device   = task.den_y;
  const auto*     dden_yx_eval_device = task.dden_yx;
  const auto*     dden_yy_eval_device = task.dden_yy;
  const auto*     dden_yz_eval_device = task.dden_yz;

  const auto*     den_x_eval_device   = task.den_x;
  const auto*     dden_xx_eval_device = task.dden_xx;
  const auto*     dden_xy_eval_device = task.dden_xy;
  const auto*     dden_xz_eval_device = task.dden_xz;

  auto*     gamma_pp_eval_device  = task.gamma_pp;
  auto*     gamma_pm_eval_device  = task.gamma_pm;
  auto*     gamma_mm_eval_device  = task.gamma_mm;

  auto*     H_z_eval_device = task.H_z;
  auto*     H_y_eval_device = task.H_y;
  auto*     H_x_eval_device = task.H_x;
  auto*     K_z_eval_device = task.K_z;
  auto*     K_y_eval_device = task.K_y;
  auto*     K_x_eval_device = task.K_x;

  const double dtolsq = 1e-24;  // TODO: make variable

  const int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if( tid < npts ) {
    const double dndz = dden_sz_eval_device[ tid ];
    const double dndy = dden_sy_eval_device[ tid ];
    const double dndx = dden_sx_eval_device[ tid ];

    const double dMzdz = dden_zz_eval_device[ tid ];
    const double dMzdy = dden_zy_eval_device[ tid ];
    const double dMzdx = dden_zx_eval_device[ tid ];

    const double dMydz = dden_yz_eval_device[ tid ];
    const double dMydy = dden_yy_eval_device[ tid ];
    const double dMydx = dden_yx_eval_device[ tid ];

    const double dMxdz = dden_xz_eval_device[ tid ];
    const double dMxdy = dden_xy_eval_device[ tid ];
    const double dMxdx = dden_xx_eval_device[ tid ];

    const auto ps = den_s_eval_device[ tid ];
    const auto pz = den_z_eval_device[ tid ];
    const auto py = den_y_eval_device[ tid ];
    const auto px = den_x_eval_device[ tid ];

    const auto mtemp = pz*pz + px*px + py*py;
    double mnorm = 0.;

    auto dels_dot_dels = dndx * dndx + dndy * dndy + dndz * dndz;
    auto delz_dot_delz = dMzdx * dMzdx + dMzdy * dMzdy + dMzdz * dMzdz;
    auto delx_dot_delx = dMxdx * dMxdx + dMxdy * dMxdy + dMxdz * dMxdz;
    auto dely_dot_dely = dMydx * dMydx + dMydy * dMydy + dMydz * dMydz;

    auto dels_dot_delz = dndx * dMzdx + dndy * dMzdy + dndz * dMzdz;
    auto dels_dot_delx = dndx * dMxdx + dndy * dMxdy + dndz * dMxdz;
    auto dels_dot_dely = dndx * dMydx + dndy * dMydy + dndz * dMydz;

    auto sum = delz_dot_delz + delx_dot_delx + dely_dot_dely;
    auto s_sum =
        dels_dot_delz * pz + dels_dot_delx * px + dels_dot_dely * py;

    auto sqsum2 =
        sqrt(dels_dot_delz * dels_dot_delz + dels_dot_delx * dels_dot_delx +
             dels_dot_dely * dels_dot_dely);

    double sign = 1.;
    if( signbit(s_sum)) 
      sign = -1.;

    if (mtemp > dtolsq) {
      mnorm = sqrt(mtemp);
      K_z_eval_device[ tid ] = pz / mnorm;
      K_y_eval_device[ tid ] = py / mnorm;
      K_x_eval_device[ tid ] = px / mnorm;
      H_z_eval_device[ tid ] = sign * dels_dot_delz / sqsum2;
      H_y_eval_device[ tid ] = sign * dels_dot_dely / sqsum2;
      H_x_eval_device[ tid ] = sign * dels_dot_delx / sqsum2;
    }
    else {
      mnorm = (1. / 3.) * (px + py + pz);
      K_z_eval_device[ tid ] = 1. / 3.;
      K_y_eval_device[ tid ] = 1. / 3.;
      K_x_eval_device[ tid ] = 1. / 3.;

      H_z_eval_device[ tid ] = sign / 3.;
      H_y_eval_device[ tid ] = sign / 3.;
      H_x_eval_device[ tid ] = sign / 3.;
    }

    gamma_pp_eval_device[ tid ] = 0.25*(dels_dot_dels + sum) + 0.5*sign*sqsum2;
    gamma_pm_eval_device[ tid ] = 0.25*(dels_dot_dels - sum);
    gamma_mm_eval_device[ tid ] = 0.25*(dels_dot_dels + sum) - 0.5*sign*sqsum2;

    den_s_eval_device[ tid ] = 0.5*(ps + mnorm);
    den_z_eval_device[ tid ] = 0.5*(ps - mnorm);

  }

}




void eval_uvars_lda_( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( npts_max,  threads.x ),
               1,
               ntasks );
  switch ( ks_scheme ) {
    case RKS:
      break;
    case UKS:
      eval_uvars_lda_uks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case GKS:
      eval_uvars_lda_gks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate UV vars" );
  }

}



void eval_uvars_gga_( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( npts_max,  threads.x ),
               1,
               ntasks );
  switch ( ks_scheme ) {
    case RKS:
      eval_uvars_gga_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case UKS:
      eval_uvars_gga_uks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case GKS:
      eval_uvars_gga_gks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate UV vars" );
  }

}









template <density_id den_select>
__global__ void eval_vvar_grad_kern( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  double* den_eval_device   = nullptr;
  double* den_x_eval_device = nullptr;
  double* den_y_eval_device = nullptr;
  double* den_z_eval_device = nullptr;

  constexpr auto warp_size = cuda::warp_size;

  if constexpr (den_select == DEN_S) {
    den_eval_device   = task.den_s;
    den_x_eval_device = task.dden_sx;
    den_y_eval_device = task.dden_sy;
    den_z_eval_device = task.dden_sz;
  }
  if constexpr (den_select == DEN_Z) {
    den_eval_device   = task.den_z;
    den_x_eval_device = task.dden_zx;
    den_y_eval_device = task.dden_zy;
    den_z_eval_device = task.dden_zz;
  }
  if constexpr (den_select == DEN_Y) {
    den_eval_device   = task.den_y;
    den_x_eval_device = task.dden_yx;
    den_y_eval_device = task.dden_yy;
    den_z_eval_device = task.dden_yz;
  }
  if constexpr (den_select == DEN_X) {
    den_eval_device   = task.den_x;
    den_x_eval_device = task.dden_xx;
    den_y_eval_device = task.dden_xy;
    den_z_eval_device = task.dden_xz;
  }

  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  const auto* den_basis_prod_device = task.zmat;
  
  __shared__ double den_shared[4][warp_size][GGA_KERNEL_SM_BLOCK_Y+1];

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
        den_shared[3][threadIdx.x][sm_y] = 0.;

        if (bid_y + threadIdx.x < npts and bid_x + sm_y < nbf) { 
          const double* db_col   = den_basis_prod_device + (bid_x + sm_y)*npts;
          const double* bf_col   = basis_eval_device     + (bid_x + sm_y)*npts;
          const double* bf_x_col = dbasis_x_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_y_col = dbasis_y_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_z_col = dbasis_z_eval_device  + (bid_x + sm_y)*npts;

          den_shared[0][threadIdx.x][sm_y] = bf_col  [ bid_y + threadIdx.x ] * db_col[ bid_y + threadIdx.x ];
          den_shared[1][threadIdx.x][sm_y] = bf_x_col[ bid_y + threadIdx.x ] * db_col[ bid_y + threadIdx.x ];
          den_shared[2][threadIdx.x][sm_y] = bf_y_col[ bid_y + threadIdx.x ] * db_col[ bid_y + threadIdx.x ];
          den_shared[3][threadIdx.x][sm_y] = bf_z_col[ bid_y + threadIdx.x ] * db_col[ bid_y + threadIdx.x ];
        }
      }
      __syncthreads();


      for (int sm_y = threadIdx.y; sm_y < GGA_KERNEL_SM_BLOCK_Y; sm_y += blockDim.y) {
        const int tid_y = bid_y + sm_y;
        register double den_reg = den_shared[0][sm_y][threadIdx.x];
        register double dx_reg  = den_shared[1][sm_y][threadIdx.x];
        register double dy_reg  = den_shared[2][sm_y][threadIdx.x];
        register double dz_reg  = den_shared[3][sm_y][threadIdx.x];

        // Warp blocks are stored col major
        den_reg =     cuda::warp_reduce_sum<warp_size>( den_reg );
        dx_reg  = 2. * cuda::warp_reduce_sum<warp_size>( dx_reg );
        dy_reg  = 2. * cuda::warp_reduce_sum<warp_size>( dy_reg );
        dz_reg  = 2. * cuda::warp_reduce_sum<warp_size>( dz_reg );


        if( threadIdx.x == 0 and tid_y < npts ) {
          atomicAdd( den_eval_device   + tid_y, den_reg );
          atomicAdd( den_x_eval_device + tid_y, dx_reg  );
          atomicAdd( den_y_eval_device + tid_y, dy_reg  );
          atomicAdd( den_z_eval_device + tid_y, dz_reg  );
        }
      }
      __syncthreads();
    }
  }

}




template <density_id den_select>
__global__ void eval_vvar_kern( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  double* den_eval_device   = nullptr;
  // use the "U" variable (+/- for UKS) even though at this point the density (S/Z) is stored
  if constexpr (den_select == DEN_S) den_eval_device = task.den_s;
  if constexpr (den_select == DEN_Z) den_eval_device = task.den_z;
  if constexpr (den_select == DEN_Y) den_eval_device = task.den_y;
  if constexpr (den_select == DEN_X) den_eval_device = task.den_x;

  const auto* basis_eval_device = task.bf;

  const auto* den_basis_prod_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  register double den_reg = 0.;

  if( tid_x < nbf and tid_y < npts ) {

    const double* bf_col   = basis_eval_device     + tid_x*npts;
    const double* db_col   = den_basis_prod_device + tid_x*npts;

    den_reg = bf_col[ tid_y ]   * db_col[ tid_y ];

  }

  // Warp blocks are stored col major
  constexpr auto warp_size = cuda::warp_size;
  //constexpr auto max_warps_per_thread_block = cuda::max_warps_per_thread_block;
  den_reg = cuda::warp_reduce_sum<warp_size>( den_reg );


  if( threadIdx.x == 0 and tid_y < npts ) {
    atomicAdd( den_eval_device   + tid_y, den_reg );
  }
  

}






void eval_vvar_( size_t ntasks, int32_t nbf_max, int32_t npts_max, bool do_grad, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
  switch( den_select ) {
    case DEN_S: 
      if (do_grad)  eval_vvar_grad_kern<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      else          eval_vvar_kern<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case DEN_Z: 
      if (do_grad)  eval_vvar_grad_kern<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      else          eval_vvar_kern<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case DEN_Y: 
      if (do_grad)  eval_vvar_grad_kern<DEN_Y><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      else          eval_vvar_kern<DEN_Y><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    case DEN_X: 
      if (do_grad)  eval_vvar_grad_kern<DEN_X><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      else          eval_vvar_kern<DEN_X><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "eval_vvar called with improper density selected" );
  }

}










}
