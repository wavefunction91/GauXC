#include "hip/hip_runtime.h"
#include "hip_eval_denvars.hpp"
#include "hip_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "hip_device_properties.hpp"

namespace GauXC      {
namespace integrator {
namespace hip       {
using namespace GauXC::hip;

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

  double den_reg = 0.;

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


  for( int  ipt = blockIdx.y * blockDim.y + threadIdx.y;
            ipt < npts;
	    ipt += blockDim.y * gridDim.y ) {

    double den = 0.;
    double dx  = 0.;
    double dy  = 0.;
    double dz  = 0.;
    
    for( int ibf_st = 0; ibf_st < nbf; ibf_st += warp_size ) {

      double den_reg = 0.;
      double dx_reg  = 0.;
      double dy_reg  = 0.;
      double dz_reg  = 0.;

      int ibf = ibf_st + threadIdx.x;
      if( ibf < nbf ) {
        const double* bf_col   = basis_eval_device     + ibf*npts;
        const double* bf_x_col = dbasis_x_eval_device  + ibf*npts;
        const double* bf_y_col = dbasis_y_eval_device  + ibf*npts;
        const double* bf_z_col = dbasis_z_eval_device  + ibf*npts;
        const double* db_col   = den_basis_prod_device + ibf*npts;

        den_reg = bf_col[ ipt ]   * db_col[ ipt ];
        dx_reg  = bf_x_col[ ipt ] * db_col[ ipt ];
        dy_reg  = bf_y_col[ ipt ] * db_col[ ipt ];
        dz_reg  = bf_z_col[ ipt ] * db_col[ ipt ];
      }

      den += 2 * warpReduceSum( den_reg );
      dx  += 4 * warpReduceSum( dx_reg );
      dy  += 4 * warpReduceSum( dy_reg );
      dz  += 4 * warpReduceSum( dz_reg );
      
    }

    if( threadIdx.x == 0 ) {
      den_eval_device   [ipt] = den;
      den_x_eval_device [ipt] = dx ;
      den_y_eval_device [ipt] = dy ;
      den_z_eval_device [ipt] = dz ;
    }
    //__sync_warp();

  }


/*
  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  double den_reg = 0.;
  double dx_reg  = 0.;
  double dy_reg  = 0.;
  double dz_reg  = 0.;

  if( tid_x < nbf and tid_y < npts ) {

    const double* bf_col   = basis_eval_device     + tid_x*npts;
    const double* bf_x_col = dbasis_x_eval_device  + tid_x*npts;
    const double* bf_y_col = dbasis_y_eval_device  + tid_x*npts;
    const double* bf_z_col = dbasis_z_eval_device  + tid_x*npts;
    const double* db_col   = den_basis_prod_device + tid_x*npts;

    den_reg = bf_col[ tid_y ]   * db_col[ tid_y ];
    dx_reg  = bf_x_col[ tid_y ] * db_col[ tid_y ];
    dy_reg  = bf_y_col[ tid_y ] * db_col[ tid_y ];
    dz_reg  = bf_z_col[ tid_y ] * db_col[ tid_y ];

  }

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
*/
  

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
                            hipStream_t     stream ) {

  dim3 threads(warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( util::div_ceil( max_nbf , threads.x ),
               util::div_ceil( max_npts , threads.y ),
               ntasks );

  hipLaunchKernelGGL(eval_uvars_lda_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}

template <typename T>
void eval_uvars_gga_device( size_t           ntasks,
                            size_t           max_nbf,
                            size_t           max_npts,
                            XCTaskDevice<T>* tasks_device,
                            hipStream_t     stream ) {

  dim3 threads(warp_size, max_warps_per_thread_block, 1);
  dim3 blocks( 1, 8, ntasks );

  hipLaunchKernelGGL(eval_uvars_gga_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}
 

template <typename T>
void eval_vvars_gga_device( size_t       npts,
                            const T*     den_x_device,
                            const T*     den_y_device,
                            const T*     den_z_device,
                                  T*     gamma_device,
                            hipStream_t stream ) {

  dim3 threads( max_threads_per_thread_block );
  dim3 blocks( util::div_ceil( npts, threads.x ) );

  hipLaunchKernelGGL(eval_vvars_gga_kernel, dim3(blocks), dim3(threads), 0, stream , 
    npts, den_x_device, den_y_device, den_z_device, gamma_device
  );

}
                          














template
void eval_uvars_lda_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            hipStream_t          stream );

template
void eval_uvars_gga_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            hipStream_t          stream );

template
void eval_vvars_gga_device( size_t            npts,
                            const double*     den_x_device,
                            const double*     den_y_device,
                            const double*     den_z_device,
                                  double*     gamma_device,
                            hipStream_t      stream );

}
}
}
