#include "cuda/cuda_eval_denvars.hpp"
#include "cuda/cuda_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "cuda/cuda_device_properties.hpp"

namespace GauXC      {
namespace integrator {
namespace cuda       {

using namespace GauXC::cuda;

template <typename T>
__global__ void eval_uvars_lda_kernel( size_t           ntasks,
                                       XCTaskDevice<T>* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.nbe;

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
  den_reg = 2 * warpReduceSum( den_reg );


  if( threadIdx.x == 0 and tid_y < npts ) {
    atomicAdd( den_eval_device   + tid_y, den_reg );
  }
  

}



#define GGA_KERNEL_SM_BLOCK_Y 32

template <typename T>
__global__ void eval_uvars_gga_kernel( size_t           ntasks,
                                       XCTaskDevice<T>* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.nbe;

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
        den_reg = 2 * warpReduceSum( den_reg );
        dx_reg  = 4 * warpReduceSum( dx_reg );
        dy_reg  = 4 * warpReduceSum( dy_reg );
        dz_reg  = 4 * warpReduceSum( dz_reg );


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


template <typename T>
__global__ void eval_vvars_gga_kernel( 
  size_t   npts,
  const T* den_x_eval_device,
  const T* den_y_eval_device,
  const T* den_z_eval_device,
        T* gamma_eval_device
) {

  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < npts ) {

    const double dx = den_x_eval_device[ tid ];
    const double dy = den_y_eval_device[ tid ];
    const double dz = den_z_eval_device[ tid ];

    gamma_eval_device[tid] = dx*dx + dy*dy + dz*dz;

  }

}


template <typename T>
void eval_uvars_lda_device( size_t           ntasks,
                            size_t           max_nbf,
                            size_t           max_npts,
                            XCTaskDevice<T>* tasks_device,
                            cudaStream_t     stream ) {

  dim3 threads(warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( max_nbf , threads.x ),
               util::div_ceil( max_npts , threads.y ),
               ntasks );

  eval_uvars_lda_kernel<<< blocks, threads, 0, stream >>>( ntasks, tasks_device );

}

template <typename T>
void eval_uvars_gga_device( size_t           ntasks,
                            size_t           max_nbf,
                            size_t           max_npts,
                            XCTaskDevice<T>* tasks_device,
                            cudaStream_t     stream ) {

  dim3 threads( warp_size, max_warps_per_thread_block / 2, 1 );
  dim3 blocks( std::min(int64_t(4), util::div_ceil( max_nbf, 4 )),
               std::min(int64_t(16), util::div_ceil( max_nbf, 16 )),
               ntasks );

  eval_uvars_gga_kernel<<< blocks, threads, 0, stream >>>( ntasks, tasks_device );

}
 

template <typename T>
void eval_vvars_gga_device( size_t       npts,
                            const T*     den_x_device,
                            const T*     den_y_device,
                            const T*     den_z_device,
                                  T*     gamma_device,
                            cudaStream_t stream ) {

  dim3 threads( max_threads_per_thread_block );
  dim3 blocks( util::div_ceil( npts, threads.x ) );

  eval_vvars_gga_kernel<<< blocks, threads, 0, stream >>>(
    npts, den_x_device, den_y_device, den_z_device, gamma_device
  );

}
                          














template
void eval_uvars_lda_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            cudaStream_t          stream );

template
void eval_uvars_gga_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            cudaStream_t          stream );

template
void eval_vvars_gga_device( size_t            npts,
                            const double*     den_x_device,
                            const double*     den_y_device,
                            const double*     den_z_device,
                                  double*     gamma_device,
                            cudaStream_t      stream );

}
}
}
