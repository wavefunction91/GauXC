#include "hip/hip_runtime.h"
#include "device/hip/hip_zmat.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device/hip/hip_device_properties.hpp"

namespace GauXC      {
namespace integrator {
namespace hip       {

using namespace GauXC::hip;


template <typename T>
__global__ void zmat_lda_kernel( size_t           ntasks,
                                 XCTaskDevice<T>* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.nbe;
  const auto* vrho_device    = task.vrho;

  const auto* basis_eval_device = task.bf;

  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;
    const double fact = 0.5 * vrho_device[tid_x];

    z_matrix_device[ ibfoff ] = fact * basis_eval_device[ ibfoff ];

  }

}




template <typename T>
void zmat_lda_hip( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    hipStream_t     stream ) {


  dim3 threads(warp_size,max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  hipLaunchKernelGGL(zmat_lda_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}

template
void zmat_lda_hip( size_t                ntasks,
                    int32_t               max_nbf,
                    int32_t               max_npts,
                    XCTaskDevice<double>* tasks_device,
                    hipStream_t          stream ); 




template <typename T>
__global__ void zmat_gga_kernel( size_t           ntasks,
                                 XCTaskDevice<T>* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.nbe;
  const auto* vrho_device    = task.vrho;
  const auto* vgamma_device  = task.vgamma;
  const auto* den_x_eval_device = task.ddenx;
  const auto* den_y_eval_device = task.ddeny;
  const auto* den_z_eval_device = task.ddenz;

  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;
    const double fact_1 = 0.5 * vrho_device[tid_x]  ;
    const double fact_2 = 2.0 * vgamma_device[tid_x];

    const double dx = den_x_eval_device[ tid_x ] * dbasis_x_eval_device[ ibfoff ];
    const double dy = den_y_eval_device[ tid_x ] * dbasis_y_eval_device[ ibfoff ];
    const double dz = den_z_eval_device[ tid_x ] * dbasis_z_eval_device[ ibfoff ];

    z_matrix_device[ ibfoff ] = 
      fact_1 * basis_eval_device[ ibfoff ] + fact_2 * ( dx + dy + dz ); 

  }
}

template <typename T>
void zmat_gga_hip( size_t           ntasks,
                    int32_t          max_nbf,
                    int32_t          max_npts,
                    XCTaskDevice<T>* tasks_device,
                    hipStream_t     stream ) {


  dim3 threads(warp_size,max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  hipLaunchKernelGGL(zmat_gga_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}
template
void zmat_gga_hip( size_t                ntasks,
                    int32_t               max_nbf,
                    int32_t               max_npts,
                    XCTaskDevice<double>* tasks_device,
                    hipStream_t          stream ); 
              
}
}
}

