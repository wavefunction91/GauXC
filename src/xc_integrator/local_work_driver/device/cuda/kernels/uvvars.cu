/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/uvvars.hpp"
#include "cuda_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "uvvars_lda.hpp"
#include "uvvars_gga.hpp"
#include "uvvars_mgga.hpp"

namespace GauXC {

#define EVAL_UVARS_KERNEL(xc_approx) \
  cudaStream_t stream = queue.queue_as<util::cuda_stream>();  \
  switch ( ks_scheme ) { \
    case RKS: \
      eval_uvars_##xc_approx##_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    case UKS: \
      eval_uvars_##xc_approx##_uks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    case GKS: \
      eval_uvars_##xc_approx##_gks_kernel<<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate U vars" ); \
  } 

#define EVAL_VVARS_KERNEL(xc_approx) \
  cudaStream_t stream = queue.queue_as<util::cuda_stream>();  \
  switch ( den_select ) { \
    case DEN_S: \
      eval_vvar_##xc_approx##_kern<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    case DEN_Z: \
      eval_vvar_##xc_approx##_kern<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    case DEN_Y: \
      eval_vvar_##xc_approx##_kern<DEN_Y><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    case DEN_X: \
      eval_vvar_##xc_approx##_kern<DEN_X><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); \
      break; \
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate V vars" ); \
  } 

void eval_uvars_lda( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  dim3 threads( cuda::max_warps_per_thread_block * cuda::warp_size, 1, 1 );
  dim3 blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks ); 
  EVAL_UVARS_KERNEL(lda);
}

void eval_vvars_lda( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
  EVAL_VVARS_KERNEL(lda);

}






void eval_uvars_gga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  dim3 threads( GGA_KERNEL_SM_WARPS * cuda::warp_size, 1, 1 );
  dim3 blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks ); 
  EVAL_UVARS_KERNEL(gga);
}

void eval_vvars_gga( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
  EVAL_VVARS_KERNEL(gga);
}


void eval_uvars_mgga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>(); 

  // Evaluate GAMMA
  eval_uvars_gga(ntasks, npts_max, ks_scheme, device_tasks, queue);

  if(ks_scheme == RKS) {
    return; // Nothing left to do
  } else if(ks_scheme == UKS) {
    dim3 threads( cuda::max_warps_per_thread_block * cuda::warp_size, 1, 1 );
    dim3 blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks ); 
    if(need_lapl) {
      eval_uvars_mgga_uks_kernel<true><<<blocks, threads, 0, stream>>>(ntasks, device_tasks);
    } else {
      eval_uvars_mgga_uks_kernel<false><<<blocks, threads, 0, stream>>>(ntasks, device_tasks);
    }
  } else {
    GAUXC_GENERIC_EXCEPTION("GKS + MGGA + DEVICE NYI");
  }

}

void eval_vvars_mgga( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {

  // First evaluate GGA variables
  eval_vvars_gga(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);

  dim3 threads( cuda::warp_size, cuda::max_warps_per_thread_block, 1 );
  dim3 blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  switch ( den_select ) {
    case DEN_S:
      if(need_lapl)
        eval_vvar_mgga_kern<DEN_S,true><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); 
      else
        eval_vvar_mgga_kern<DEN_S,false><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); 
      break;
    case DEN_Z:
      if(need_lapl)
        eval_vvar_mgga_kern<DEN_Z,true><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); 
      else
        eval_vvar_mgga_kern<DEN_Z,false><<< blocks, threads, 0, stream >>>( ntasks, device_tasks ); 
      break;
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate V vars" ); \
  } 
}

}
