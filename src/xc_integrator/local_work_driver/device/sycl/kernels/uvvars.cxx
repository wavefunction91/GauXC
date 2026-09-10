/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/uvvars.hpp"
#include "device_specific/sycl_util.hpp"
#include "sycl_launch.hpp"
#include <gauxc/util/div_ceil.hpp>

#include "uvvars_lda.hpp"
#include "uvvars_gga.hpp"
#include "uvvars_mgga.hpp"

namespace GauXC {

#define EVAL_UVARS_KERNEL(xc_approx) \
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>(); \
  switch ( ks_scheme ) { \
    case RKS: \
      sycl::launch_kernel( stream, blocks, threads, "eval_uvars_" #xc_approx "_rks_kernel", \
        [=](){ eval_uvars_##xc_approx##_rks_kernel( ntasks, device_tasks ); } ); \
      break; \
    case UKS: \
      sycl::launch_kernel( stream, blocks, threads, "eval_uvars_" #xc_approx "_uks_kernel", \
        [=](){ eval_uvars_##xc_approx##_uks_kernel( ntasks, device_tasks ); } ); \
      break; \
    case GKS: \
      sycl::launch_kernel( stream, blocks, threads, "eval_uvars_" #xc_approx "_gks_kernel", \
        [=](){ eval_uvars_##xc_approx##_gks_kernel( ntasks, device_tasks ); } ); \
      break; \
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate U vars" ); \
  }


#define EVAL_TMAT_KERNEL(xc_approx) \
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>(); \
  switch ( ks_scheme ) { \
    case RKS: \
      sycl::launch_kernel( stream, blocks, threads, "eval_tmat_" #xc_approx "_rks_kernel", \
        [=](){ eval_tmat_##xc_approx##_rks_kernel( ntasks, device_tasks ); } ); \
      break; \
    case UKS: \
      sycl::launch_kernel( stream, blocks, threads, "eval_tmat_" #xc_approx "_uks_kernel", \
        [=](){ eval_tmat_##xc_approx##_uks_kernel( ntasks, device_tasks ); } ); \
      break; \
    case GKS: \
      GAUXC_GENERIC_EXCEPTION( "GKS + evaluate trial U vars NYI" ); \
      break; \
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate U vars" ); \
  }


#define EVAL_VVARS_KERNEL(xc_approx) \
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>(); \
  switch ( den_select ) { \
    case DEN_S: \
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_" #xc_approx "_kern", \
        [=](){ eval_vvar_##xc_approx##_kern<trial,DEN_S>( ntasks, device_tasks ); } ); \
      break; \
    case DEN_Z: \
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_" #xc_approx "_kern", \
        [=](){ eval_vvar_##xc_approx##_kern<trial,DEN_Z>( ntasks, device_tasks ); } ); \
      break; \
    case DEN_Y: \
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_" #xc_approx "_kern", \
        [=](){ eval_vvar_##xc_approx##_kern<trial,DEN_Y>( ntasks, device_tasks ); } ); \
      break; \
    case DEN_X: \
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_" #xc_approx "_kern", \
        [=](){ eval_vvar_##xc_approx##_kern<trial,DEN_X>( ntasks, device_tasks ); } ); \
      break; \
    default: \
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate V vars" ); \
  }

// Internal implementation with trial parameter
void eval_tmat_lda( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  sycl::launch_dim threads( sycl::max_warps_per_thread_block * sycl::warp_size, 1, 1 );
  sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );
  EVAL_TMAT_KERNEL(lda);
}

void eval_uvars_lda( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  sycl::launch_dim threads( sycl::max_warps_per_thread_block * sycl::warp_size, 1, 1 );
  sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );
  EVAL_UVARS_KERNEL(lda);
}

// Internal implementation with trial as template parameter
template<bool trial>
void eval_vvars_lda_impl( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  sycl::launch_dim threads( sycl::warp_size, sycl::max_warps_per_thread_block, 1 );
  sycl::launch_dim blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );
  EVAL_VVARS_KERNEL(lda);
}
void eval_vvars_lda( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_lda_impl<false>(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);
}
void eval_vvars_lda_trial( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_lda_impl<true>(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);
}

// Internal implementation with trial parameter
void eval_tmat_gga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  sycl::launch_dim threads( GGA_KERNEL_SM_WARPS * sycl::warp_size, 1, 1 );
  sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );
  EVAL_TMAT_KERNEL(gga);
}
void eval_uvars_gga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  XCDeviceTask* device_tasks, device_queue queue ) {
  sycl::launch_dim threads( GGA_KERNEL_SM_WARPS * sycl::warp_size, 1, 1 );
  sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );
  EVAL_UVARS_KERNEL(gga);
}


// Internal implementation with trial as template parameter
template<bool trial>
void eval_vvars_gga_impl( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();
  sycl::launch_dim threads( sycl::warp_size, sycl::max_warps_per_thread_block, 1 );
  sycl::launch_dim blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );

  switch ( den_select ) {
    case DEN_S:
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_gga_kern",
        [=](){ eval_vvar_gga_kern<trial,DEN_S>( ntasks, device_tasks ); } );
      break;
    case DEN_Z:
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_gga_kern",
        [=](){ eval_vvar_gga_kern<trial,DEN_Z>( ntasks, device_tasks ); } );
      break;
    case DEN_Y:
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_gga_kern",
        [=](){ eval_vvar_gga_kern<trial,DEN_Y>( ntasks, device_tasks ); } );
      break;
    case DEN_X:
      sycl::launch_kernel( stream, blocks, threads, "eval_vvar_gga_kern",
        [=](){ eval_vvar_gga_kern<trial,DEN_X>( ntasks, device_tasks ); } );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate V vars" );
  }
}
void eval_vvars_gga( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_gga_impl<false>(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);
}
void eval_vvars_gga_trial( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_gga_impl<true>(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);
}

// Internal implementation with trial parameter
void eval_tmat_mgga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {
  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  sycl::launch_dim threads( GGA_KERNEL_SM_WARPS * sycl::warp_size, 1, 1 );
  sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );

  if(need_lapl) {
    GAUXC_GENERIC_EXCEPTION("MGGA + LAPL + eval tmat NYI");
  }
  if(ks_scheme == RKS) {
      sycl::launch_kernel( stream, blocks, threads, "eval_tmat_mgga_rks_kernel",
        [=](){ eval_tmat_mgga_rks_kernel( ntasks, device_tasks ); } );
  } else if(ks_scheme == UKS) {
      sycl::launch_kernel( stream, blocks, threads, "eval_tmat_mgga_uks_kernel",
        [=](){ eval_tmat_mgga_uks_kernel( ntasks, device_tasks ); } );
  } else {
    GAUXC_GENERIC_EXCEPTION("GKS + MGGA + DEVICE NYI");
  }
}

void eval_uvars_mgga( size_t ntasks, int32_t npts_max, integrator_ks_scheme ks_scheme,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();

  // Evaluate GAMMA
  eval_uvars_gga(ntasks, npts_max, ks_scheme, device_tasks, queue);

  if(ks_scheme == RKS) {
    return; // Nothing left to do
  } else if(ks_scheme == UKS) {
    sycl::launch_dim threads( sycl::max_warps_per_thread_block * sycl::warp_size, 1, 1 );
    sycl::launch_dim blocks( util::div_ceil( npts_max,  threads.x ), 1, ntasks );
    if(need_lapl) {
      sycl::launch_kernel( stream, blocks, threads, "eval_uvars_mgga_uks_kernel",
        [=](){ eval_uvars_mgga_uks_kernel<true>( ntasks, device_tasks ); } );
    } else {
      sycl::launch_kernel( stream, blocks, threads, "eval_uvars_mgga_uks_kernel",
        [=](){ eval_uvars_mgga_uks_kernel<false>( ntasks, device_tasks ); } );
    }
  } else {
    GAUXC_GENERIC_EXCEPTION("GKS + MGGA + DEVICE NYI");
  }

}

// Internal implementation with trial as template parameter
template<bool trial>
void eval_vvars_mgga_impl( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {
  // First evaluate GGA variables
  eval_vvars_gga_impl<trial>(ntasks, nbf_max, npts_max, den_select, device_tasks, queue);

  ::sycl::queue& stream = queue.queue_as<util::sycl_queue>();
  sycl::launch_dim threads( sycl::warp_size, sycl::max_warps_per_thread_block, 1 );
  sycl::launch_dim blocks( util::div_ceil( nbf_max,  threads.x ),
               util::div_ceil( npts_max, threads.y ),
               ntasks );

  switch ( den_select ) {
    case DEN_S:
      if (need_lapl) {
        sycl::launch_kernel( stream, blocks, threads, "eval_vvar_mgga_kern",
          [=](){ eval_vvar_mgga_kern<trial,DEN_S,true>( ntasks, device_tasks ); } );
      } else {
        sycl::launch_kernel( stream, blocks, threads, "eval_vvar_mgga_kern",
          [=](){ eval_vvar_mgga_kern<trial,DEN_S,false>( ntasks, device_tasks ); } );
      }
      break;
    case DEN_Z:
      if (need_lapl) {
        sycl::launch_kernel( stream, blocks, threads, "eval_vvar_mgga_kern",
          [=](){ eval_vvar_mgga_kern<trial,DEN_Z,true>( ntasks, device_tasks ); } );
      } else {
        sycl::launch_kernel( stream, blocks, threads, "eval_vvar_mgga_kern",
          [=](){ eval_vvar_mgga_kern<trial,DEN_Z,false>( ntasks, device_tasks ); } );
      }
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "Unexpected KS scheme when attempting to evaluate V vars" );
  }
}
void eval_vvars_mgga( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_mgga_impl<false>(ntasks, nbf_max, npts_max, den_select, need_lapl, device_tasks, queue);
}
void eval_vvars_mgga_trial( size_t ntasks, int32_t nbf_max, int32_t npts_max, density_id den_select,
  bool need_lapl, XCDeviceTask* device_tasks, device_queue queue ) {
  eval_vvars_mgga_impl<true>(ntasks, nbf_max, npts_max, den_select, need_lapl, device_tasks, queue);
}

}
