/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "device/common/zmat_vxc.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"
#include "device_specific/cuda_device_constants.hpp"

namespace GauXC {


__global__ void zmat_lda_vxc_rks_kernel( size_t        ntasks,
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













template<density_id den_selector>
__global__ void zmat_lda_vxc_uks_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;


  const auto* basis_eval_device = task.bf;


  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;
    const double factp = 0.5 * vrho_pos_device[tid_x];
    const double factm = 0.5 * vrho_neg_device[tid_x];
    if constexpr ( den_selector == DEN_S ) // positive density
      z_matrix_device[ ibfoff ] = 0.5*(factp * basis_eval_device[ ibfoff ] + factm * basis_eval_device[ ibfoff ]);
    if constexpr ( den_selector == DEN_Z ) // negative density
      z_matrix_device[ ibfoff ] = 0.5*(factp * basis_eval_device[ ibfoff ] - factm * basis_eval_device[ ibfoff ]);
  }

}

template<density_id den_selector>
__global__ void zmat_lda_vxc_gks_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;
  const double* KZ_device          = task.K_z;
  const double* KY_device          = task.K_y;
  const double* KX_device          = task.K_x;


  const auto* basis_eval_device = task.bf;


  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;
    const double factp = 0.5 * vrho_pos_device[tid_x];
    const double factm = 0.5 * vrho_neg_device[tid_x];
    const double factk = 0.5 * (factp - factm);

    if constexpr ( den_selector == DEN_S )
      z_matrix_device[ ibfoff ] = 0.5*(factp * basis_eval_device[ ibfoff ] + factm * basis_eval_device[ ibfoff ]);
    if constexpr ( den_selector == DEN_Z )
      z_matrix_device[ ibfoff ] = KZ_device[ ibfoff ] * factk * basis_eval_device[ ibfoff ];
    if constexpr ( den_selector == DEN_Y )
      z_matrix_device[ ibfoff ] = KY_device[ ibfoff ] * factk * basis_eval_device[ ibfoff ];
    if constexpr ( den_selector == DEN_X )
      z_matrix_device[ ibfoff ] = KX_device[ ibfoff ] * factk * basis_eval_device[ ibfoff ];
  }

}









__global__ void zmat_gga_vxc_rks_kernel( size_t        ntasks,
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








template<density_id den_selector>
__global__ void zmat_gga_vxc_uks_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;
  const double* vgamma_pp_device   = task.vgamma_pp;
  const double* vgamma_pm_device   = task.vgamma_pm;
  const double* vgamma_mm_device   = task.vgamma_mm;

  const auto* den_pos_x_eval_device = task.dden_sx;
  const auto* den_pos_y_eval_device = task.dden_sy;
  const auto* den_pos_z_eval_device = task.dden_sz;
  const auto* den_neg_x_eval_device = task.dden_zx;
  const auto* den_neg_y_eval_device = task.dden_zy;
  const auto* den_neg_z_eval_device = task.dden_zz;


  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;

    const double factp = 0.25 * vrho_pos_device[tid_x];
    const double factm = 0.25 * vrho_neg_device[tid_x];
    
    const auto gga_fact_pp  = vgamma_pp_device[tid_x];
    const auto gga_fact_pm  = vgamma_pm_device[tid_x];
    const auto gga_fact_mm  = vgamma_mm_device[tid_x];
    
    const auto gga_fact_1 = 0.5*(gga_fact_pp + gga_fact_pm + gga_fact_mm);
    const auto gga_fact_2 = 0.5*(gga_fact_pp - gga_fact_mm);
    const auto gga_fact_3 = 0.5*(gga_fact_pp - gga_fact_pm + gga_fact_mm);

    if constexpr ( den_selector == DEN_S ) {
      const auto x_fact = gga_fact_1 * den_pos_x_eval_device[ tid_x ] + gga_fact_2 * den_neg_x_eval_device[ tid_x ];
      const auto y_fact = gga_fact_1 * den_pos_y_eval_device[ tid_x ] + gga_fact_2 * den_neg_y_eval_device[ tid_x ];
      const auto z_fact = gga_fact_1 * den_pos_z_eval_device[ tid_x ] + gga_fact_2 * den_neg_z_eval_device[ tid_x ];
      
      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                  + (factp + factm) * basis_eval_device[ ibfoff ];

    }
    if constexpr ( den_selector == DEN_Z ) {
      const auto x_fact = gga_fact_3 * den_neg_x_eval_device[ tid_x ] + gga_fact_2 * den_pos_x_eval_device[ tid_x ];
      const auto y_fact = gga_fact_3 * den_neg_y_eval_device[ tid_x ] + gga_fact_2 * den_pos_y_eval_device[ tid_x ];
      const auto z_fact = gga_fact_3 * den_neg_z_eval_device[ tid_x ] + gga_fact_2 * den_pos_z_eval_device[ tid_x ];

      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ] 
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ]
                                  + (factp - factm) * basis_eval_device[ ibfoff ];
    }
  }
}




template<density_id den_selector>
__global__ void zmat_gga_vxc_gks_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;
  const double* vgamma_pp_device   = task.vgamma_pp;
  const double* vgamma_pm_device   = task.vgamma_pm;
  const double* vgamma_mm_device   = task.vgamma_mm;

  const double* Kz_device          = task.K_z;
  const double* Ky_device          = task.K_y;
  const double* Kx_device          = task.K_x;
  const double* Hz_device          = task.H_z;
  const double* Hy_device          = task.H_y;
  const double* Hx_device          = task.H_x;

  const auto* dden_sx_eval_device = task.dden_sx;
  const auto* dden_sy_eval_device = task.dden_sy;
  const auto* dden_sz_eval_device = task.dden_sz;
  const auto* dden_zx_eval_device = task.dden_zx;
  const auto* dden_zy_eval_device = task.dden_zy;
  const auto* dden_zz_eval_device = task.dden_zz;
  const auto* dden_yx_eval_device = task.dden_yx;
  const auto* dden_yy_eval_device = task.dden_yy;
  const auto* dden_yz_eval_device = task.dden_yz;
  const auto* dden_xx_eval_device = task.dden_xx;
  const auto* dden_xy_eval_device = task.dden_xy;
  const auto* dden_xz_eval_device = task.dden_xz;


  const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;

  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;

    const double fact_p =  0.5*vrho_pos_device[tid_x];
    const double fact_m =  0.5*vrho_neg_device[tid_x];
    
    const auto gga_fact_pp  = vgamma_pp_device[tid_x];
    const auto gga_fact_pm  = vgamma_pm_device[tid_x];
    const auto gga_fact_mm  = vgamma_mm_device[tid_x];
    
    const auto gga_fact_1 = 0.5*(gga_fact_pp + gga_fact_pm + gga_fact_mm);
    const auto gga_fact_2 = 0.5*(gga_fact_pp - gga_fact_mm);
    const auto gga_fact_3 = 0.5*(gga_fact_pp - gga_fact_pm + gga_fact_mm);

    if constexpr ( den_selector == DEN_S ) {
      const auto s_fact = 0.5 * (fact_p + fact_m);
      const auto x_fact = gga_fact_1 * dden_sx_eval_device[ tid_x ]
                        + gga_fact_2 * (Hz_device[ tid_x ] * dden_zx_eval_device[ tid_x ]
                                     +  Hy_device[ tid_x ] * dden_yx_eval_device[ tid_x ]
                                     +  Hx_device[ tid_x ] * dden_xx_eval_device[ tid_x ] );
      const auto y_fact = gga_fact_1 * dden_sy_eval_device[ tid_x ]
                        + gga_fact_2 * (Hz_device[ tid_x ] * dden_zy_eval_device[ tid_x ]
                                     +  Hy_device[ tid_x ] * dden_yy_eval_device[ tid_x ]
                                     +  Hx_device[ tid_x ] * dden_xy_eval_device[ tid_x ] );
      const auto z_fact = gga_fact_1 * dden_sz_eval_device[ tid_x ]
                        + gga_fact_2 * (Hz_device[ tid_x ] * dden_zz_eval_device[ tid_x ]
                                     +  Hy_device[ tid_x ] * dden_yz_eval_device[ tid_x ]
                                     +  Hx_device[ tid_x ] * dden_xz_eval_device[ tid_x ] );
      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                  + s_fact * basis_eval_device[ ibfoff ];

    }

    if constexpr ( den_selector == DEN_Z ) {
      const auto s_fact  = Kz_device[ tid_x ] * 0.5 * (fact_p - fact_m);
      const auto x_fact  = gga_fact_3 * dden_zx_eval_device[ tid_x ]
                        +  gga_fact_2 * Hz_device[ tid_x ] * dden_sx_eval_device[ tid_x ];
      const auto y_fact  = gga_fact_3 * dden_zy_eval_device[ tid_x ]
                        +  gga_fact_2 * Hz_device[ tid_x ] * dden_sy_eval_device[ tid_x ];
      const auto z_fact  = gga_fact_3 * dden_zz_eval_device[ tid_x ]
                        +  gga_fact_2 * Hz_device[ tid_x ] * dden_sz_eval_device[ tid_x ];
      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                  + s_fact *  basis_eval_device[ ibfoff ];
    }

    if constexpr ( den_selector == DEN_Y ) {
      const auto s_fact  = Ky_device[ tid_x ] * 0.5 * (fact_p - fact_m);
      const auto x_fact  = gga_fact_3 * dden_yx_eval_device[ tid_x ]
                        +  gga_fact_2 * Hy_device[ tid_x ] * dden_sx_eval_device[ tid_x ];
      const auto y_fact  = gga_fact_3 * dden_yy_eval_device[ tid_x ]
                        +  gga_fact_2 * Hy_device[ tid_x ] * dden_sy_eval_device[ tid_x ];
      const auto z_fact  = gga_fact_3 * dden_yz_eval_device[ tid_x ]
                        +  gga_fact_2 * Hy_device[ tid_x ] * dden_sz_eval_device[ tid_x ];
      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                  + s_fact *  basis_eval_device[ ibfoff ];
    }

    if constexpr ( den_selector == DEN_X ) {
      const auto s_fact  = Kx_device[ tid_x ] * 0.5 * (fact_p - fact_m);
      const auto x_fact  = gga_fact_3 * dden_xx_eval_device[ tid_x ]
                        +  gga_fact_2 * Hx_device[ tid_x ] * dden_sx_eval_device[ tid_x ];
      const auto y_fact  = gga_fact_3 * dden_xy_eval_device[ tid_x ]
                        +  gga_fact_2 * Hx_device[ tid_x ] * dden_sy_eval_device[ tid_x ];
      const auto z_fact  = gga_fact_3 * dden_xz_eval_device[ tid_x ]
                        +  gga_fact_2 * Hx_device[ tid_x ] * dden_sz_eval_device[ tid_x ];
      z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                  + y_fact * dbasis_y_eval_device[ ibfoff ]
                                  + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                  + s_fact *  basis_eval_device[ ibfoff ];
    }

  }
}








void zmat_lda_vxc( size_t            ntasks,
                   int32_t           max_nbf,
                   int32_t           max_npts,
                   XCDeviceTask*     tasks_device,
                   integrator_ks_scheme scheme,
                   density_id sel,
                   device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;


  dim3 threads(cuda::warp_size,cuda::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  switch( scheme ) {
    case RKS:
      zmat_lda_vxc_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      break;
    case UKS:
      if ( sel == DEN_S )       zmat_lda_vxc_uks_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Z )  zmat_lda_vxc_uks_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else GAUXC_GENERIC_EXCEPTION( "zmat_lda_vxc invalid density" );
      break;
    case GKS:
      if ( sel == DEN_S )       zmat_lda_vxc_gks_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Z )  zmat_lda_vxc_gks_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Y )  zmat_lda_vxc_gks_kernel<DEN_Y><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_X )  zmat_lda_vxc_gks_kernel<DEN_X><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else GAUXC_GENERIC_EXCEPTION( "zmat_lda_vxc invalid density" );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "zmat_lda_vxc invalid KS scheme" );
  }
}







void zmat_gga_vxc( size_t            ntasks,
                   int32_t           max_nbf,
                   int32_t           max_npts,
                   XCDeviceTask*     tasks_device,
                   integrator_ks_scheme scheme,
                   density_id sel,
                   device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;


  dim3 threads(cuda::warp_size,cuda::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
               util::div_ceil( max_nbf,  threads.y ),
               ntasks );

  switch( scheme ) {
    case RKS:
      zmat_gga_vxc_rks_kernel<<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      break;
    case UKS:
      if ( sel == DEN_S )       zmat_gga_vxc_uks_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Z )  zmat_gga_vxc_uks_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else GAUXC_GENERIC_EXCEPTION( "zmat_gga_vxc invalid density" );
      break;
    case GKS:
      if ( sel == DEN_S )       zmat_gga_vxc_gks_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Z )  zmat_gga_vxc_gks_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_Y )  zmat_gga_vxc_gks_kernel<DEN_Y><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else if ( sel == DEN_X )  zmat_gga_vxc_gks_kernel<DEN_X><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
      else GAUXC_GENERIC_EXCEPTION( "zmat_gga_vxc invalid density" );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION( "zmat_gga_vxc invalid KS scheme" );
  }
}



}

