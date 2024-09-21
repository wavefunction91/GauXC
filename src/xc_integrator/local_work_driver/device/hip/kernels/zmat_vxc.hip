/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "hip/hip_runtime.h"
#include "device/common/zmat_vxc.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/hip_util.hpp"
#include "device_specific/hip_device_constants.hpp"

namespace GauXC {


__global__ void zmat_lda_vxc_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
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




void zmat_lda_vxc( size_t            ntasks,
                   int32_t           max_nbf,
                   int32_t           max_npts,
                   XCDeviceTask*     tasks_device,
                   device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>() ;


  dim3 threads(hip::warp_size,hip::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  hipLaunchKernelGGL(zmat_lda_vxc_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}




















__global__ void zmat_gga_vxc_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const auto* vrho_device    = task.vrho;
  const auto* vgamma_device  = task.vgamma;
  const auto* den_x_eval_device = task.dden_sx;
  const auto* den_y_eval_device = task.dden_sy;
  const auto* den_z_eval_device = task.dden_sz;

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

void zmat_gga_vxc( size_t            ntasks,
                   int32_t           max_nbf,
                   int32_t           max_npts,
                   XCDeviceTask*     tasks_device,
                   device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>() ;


  dim3 threads(hip::warp_size,hip::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  hipLaunchKernelGGL(zmat_gga_vxc_kernel, dim3(blocks), dim3(threads), 0, stream ,  ntasks, tasks_device );

}
              

template <bool need_lapl>
__global__ void zmat_mgga_vxc_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const auto* vrho_device    = task.vrho;
  const auto* vgamma_device  = task.vgamma;
  const double* vlapl_device = need_lapl ? task.vlapl : nullptr;
  const auto* den_x_eval_device = task.dden_sx;
  const auto* den_y_eval_device = task.dden_sy;
  const auto* den_z_eval_device = task.dden_sz;

  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;
  const double* d2basis_lapl_eval_device = 
    need_lapl ? task.d2bflapl : nullptr;
  

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

    double val = 
      fact_1 * basis_eval_device[ ibfoff ] + fact_2 * ( dx + dy + dz ); 

    if constexpr (need_lapl) {
      val += vlapl_device[tid_x] * d2basis_lapl_eval_device[ibfoff];
    }

    z_matrix_device[ ibfoff ] = val;
  }
}

void zmat_mgga_vxc( size_t            ntasks,
                    int32_t           max_nbf,
                    int32_t           max_npts,
                    XCDeviceTask*     tasks_device,
                    bool              do_lapl,
                    device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>() ;


  dim3 threads(hip::warp_size,hip::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  if(do_lapl)
    zmat_mgga_vxc_kernel<true><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
  else
    zmat_mgga_vxc_kernel<false><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );

}


template <bool need_lapl>
__global__ void mmat_mgga_vxc_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const auto* vtau_device    = task.vtau;
  const double* vlapl_device = need_lapl ? task.vlapl : nullptr;

  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  auto* mmat_x = task.xmat_x;
  auto* mmat_y = task.xmat_y;
  auto* mmat_z = task.xmat_z;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;
    const double fact_1 = 0.25 * vtau_device[tid_x] + 
      (need_lapl ? vlapl_device[tid_x] : 0.0);

    mmat_x[ ibfoff ] = fact_1 * dbasis_x_eval_device[ ibfoff ]; 
    mmat_y[ ibfoff ] = fact_1 * dbasis_y_eval_device[ ibfoff ]; 
    mmat_z[ ibfoff ] = fact_1 * dbasis_z_eval_device[ ibfoff ]; 
  }
}

void mmat_mgga_vxc( size_t            ntasks,
                    int32_t           max_nbf,
                    int32_t           max_npts,
                    XCDeviceTask*     tasks_device,
                    bool              do_lapl,
                    device_queue queue ) {

  hipStream_t stream = queue.queue_as<util::hip_stream>() ;


  dim3 threads(hip::warp_size,hip::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  if(do_lapl)
    mmat_mgga_vxc_kernel<true><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
  else
    mmat_mgga_vxc_kernel<false><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );

}


}

