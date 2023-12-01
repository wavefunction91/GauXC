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

__global__ void eval_uvars_lda_rks_kernel( size_t        ntasks,
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



#define GGA_KERNEL_SM_BLOCK_Y 32

__global__ void eval_uvars_gga_rks_kernel( size_t           ntasks,
                                       XCDeviceTask* tasks_device ) {

  constexpr auto warp_size = cuda::warp_size;
  //constexpr auto max_warps_per_thread_block = cuda::max_warps_per_thread_block;

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  auto* den_eval_device   = task.den;
  auto* den_x_eval_device = task.ddenx;
  auto* den_y_eval_device = task.ddeny;
  auto* den_z_eval_device = task.ddenz;

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
        dx_reg  = 2 * cuda::warp_reduce_sum<warp_size>( dx_reg );
        dy_reg  = 2 * cuda::warp_reduce_sum<warp_size>( dy_reg );
        dz_reg  = 2 * cuda::warp_reduce_sum<warp_size>( dz_reg );


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

__global__ void eval_uvars_lda_uks_kernel( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  auto* den_pos_eval_device   = task.den_pos;
  auto* den_neg_eval_device   = task.den_neg;


  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;


  if( tid_y < npts ) {
		const auto ps = den_pos_eval_device[ tid_y ];
		const auto pz = den_neg_eval_device[ tid_y ];
		den_pos_eval_device[ tid_y ] = ps + pz;
		den_neg_eval_device[ tid_y ] = ps - pz;

  }
}

__global__ void eval_vvars_gga_rks_kernel( 
  size_t   npts,
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






void eval_uvvars_lda( size_t ntasks, int32_t nbf_max, int32_t npts_max, integrator_term_tracker enabled_terms,
  XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
	switch ( enabled_terms.ks_scheme ) {
		case RKS:
  		eval_uvars_lda_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
			break;
		case UKS:
  		eval_uvars_lda_uks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
			break;
		default:
			GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate UV vars" );
	}
			

}



void eval_uvvars_gga( size_t ntasks, size_t npts_total, int32_t nbf_max, 
  int32_t npts_max, integrator_term_tracker enabled_terms, XCDeviceTask* device_tasks, const double* denx, 
  const double* deny, const double* denz, double* gamma, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();

  // U Variables
  {
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block / 2, 1 );
  dim3 blocks( std::min(uint64_t(4), util::div_ceil( nbf_max, 4 )),
               std::min(uint64_t(16), util::div_ceil( nbf_max, 16 )),
               ntasks );
  eval_uvars_gga_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
  }

  // V variables
  dim3 threads( cuda::max_threads_per_thread_block );
  dim3 blocks( util::div_ceil( npts_total, threads.x ) );
  eval_vvars_gga_rks_kernel<<< blocks, threads, 0, stream >>>(
    npts_total, denx, deny, denz, gamma
  );
}










template <density_id den_select>
__global__ void eval_den_kern( size_t        ntasks,
                                       XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  double* den_eval_device   = nullptr;
  // use the "U" variable (+/- for UKS) even though at this point the density (S/Z) is stored
	if constexpr (den_select == DEN_S) den_eval_device = task.den_pos;
	if constexpr (den_select == DEN_Z) den_eval_device = task.den_neg;

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










void eval_u_den( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
	switch( den_select ) {
		case DEN_S:	
      eval_den_kern<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
			break;
		case DEN_Z:	
      eval_den_kern<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, device_tasks );
			break;
		default:
			GAUXC_GENERIC_EXCEPTION( "eval_den called with improper density selected" );
  }

}










}
