#include "device/common/onedft.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"
#include "device_specific/cuda_device_constants.hpp"

namespace GauXC {

__global__ void sz_to_ab( size_t size,
                          const double* array1, 
                          const double* array2,
                          double* result1,
                          double* result2 ) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < size) {
    double s = array1[tid];
    double z = array2[tid];
    result1[tid] = 0.5 * (s + z);
    result2[tid] = 0.5 * (s - z);
  }
}


template<density_id den_selector>
__global__ void zmat_lda_vxc_onedft_kernel( size_t        ntasks,
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
    double sign = 1.0;
    if constexpr ( den_selector == DEN_Z )  sign = -1.0;
    
    z_matrix_device[ ibfoff ] = 0.5*(factp * basis_eval_device[ ibfoff ] + sign * factm * basis_eval_device[ ibfoff ]);
  }

}

template<density_id den_selector>
__global__ void zmat_gga_vxc_onedft_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;

  const double* dden_x_grad_a   = task.gamma_pp;
  const double* dden_x_grad_b   = task.vgamma_pp;
  const double* dden_y_grad_a   = task.gamma_pm;
  const double* dden_y_grad_b   = task.vgamma_pm;
  const double* dden_z_grad_a   = task.gamma_mm;
  const double* dden_z_grad_b   = task.vgamma_mm;

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

    double sign = 1.0;

    double x_fact, y_fact, z_fact;

    if constexpr ( den_selector == DEN_S ) {
      x_fact = 0.5 * (dden_x_grad_a[tid_x] + dden_x_grad_b[tid_x]);
      y_fact = 0.5 * (dden_y_grad_a[tid_x] + dden_y_grad_b[tid_x]);
      z_fact = 0.5 * (dden_z_grad_a[tid_x] + dden_z_grad_b[tid_x]);
   }
   if constexpr ( den_selector == DEN_Z ) {
      sign = -1.0;
      x_fact = 0.5 * (dden_x_grad_a[tid_x] - dden_x_grad_b[tid_x]);
      y_fact = 0.5 * (dden_y_grad_a[tid_x] - dden_y_grad_b[tid_x]);
      z_fact = 0.5 * (dden_z_grad_a[tid_x] - dden_z_grad_b[tid_x]);
   }

    z_matrix_device[ ibfoff ] =   x_fact * dbasis_x_eval_device[ ibfoff ]      
                                + y_fact * dbasis_y_eval_device[ ibfoff ]
                                + z_fact * dbasis_z_eval_device[ ibfoff ] 
                                + (factp + sign * factm) * basis_eval_device[ ibfoff ];
  }
}

template<bool need_lapl, density_id den_selector>
__global__ void zmat_mgga_vxc_onedft_kernel( size_t        ntasks,
                                     XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;

  const double* vrho_pos_device    = task.vrho_pos;
  const double* vrho_neg_device    = task.vrho_neg;
  const double* vlapl_pos_device    = task.vlapl_pos;
  const double* vlapl_neg_device    = task.vlapl_neg;

  const double* dden_x_grad_a   = task.gamma_pp;
  const double* dden_x_grad_b   = task.vgamma_pp;
  const double* dden_y_grad_a   = task.gamma_pm;
  const double* dden_y_grad_b   = task.vgamma_pm;
  const double* dden_z_grad_a   = task.gamma_mm;
  const double* dden_z_grad_b   = task.vgamma_mm;

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
  const auto* d2basis_lapl_eval_device = task.d2bflapl;

  auto* z_matrix_device = task.zmat;

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < npts and tid_y < nbf ) {

    const size_t ibfoff = tid_y * npts + tid_x;

    const double factp = 0.25 * vrho_pos_device[tid_x];
    const double factm = 0.25 * vrho_neg_device[tid_x];
    
    double sign = 1.0;

    double x_fact, y_fact, z_fact;

    if constexpr ( den_selector == DEN_S ) {
       x_fact = 0.5 * (dden_x_grad_a[tid_x] + dden_x_grad_b[tid_x]);
       y_fact = 0.5 * (dden_y_grad_a[tid_x] + dden_y_grad_b[tid_x]);
       z_fact = 0.5 * (dden_z_grad_a[tid_x] + dden_z_grad_b[tid_x]);
    }
    if constexpr ( den_selector == DEN_Z ) {
       sign = -1.0;
       x_fact = 0.5 * (dden_x_grad_a[tid_x] - dden_x_grad_b[tid_x]);
       y_fact = 0.5 * (dden_y_grad_a[tid_x] - dden_y_grad_b[tid_x]);
       z_fact = 0.5 * (dden_z_grad_a[tid_x] - dden_z_grad_b[tid_x]);
    }

    auto val = x_fact * dbasis_x_eval_device[ ibfoff ]      
             + y_fact * dbasis_y_eval_device[ ibfoff ]
             + z_fact * dbasis_z_eval_device[ ibfoff ] 
             + (factp + sign * factm) * basis_eval_device[ ibfoff ];

    if constexpr (need_lapl) {
      const double lfactp = vlapl_pos_device[tid_x];
      const double lfactm = vlapl_neg_device[tid_x];

      val += 0.5 * (lfactp + sign * lfactm) * d2basis_lapl_eval_device[ ibfoff ];
    }

    z_matrix_device[ ibfoff ] = val;
  }
}



void zmat_onedft_vxc( size_t            ntasks,
                      int32_t           max_nbf,
                      int32_t           max_npts,
                      XCDeviceTask*     tasks_device,
                      integrator_xc_approx scheme,
                      density_id sel,
                      device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;


  dim3 threads(cuda::warp_size,cuda::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( max_npts, threads.x ),
                util::div_ceil( max_nbf,  threads.y ),
                ntasks );
  if(scheme == LDA) {
    switch(sel) {
      case DEN_S:
        zmat_lda_vxc_onedft_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
      case DEN_Z:
        zmat_lda_vxc_onedft_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
    }
  } else if(scheme == GGA) {
    switch(sel) {
      case DEN_S:
        zmat_gga_vxc_onedft_kernel<DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
      case DEN_Z:
        zmat_gga_vxc_onedft_kernel<DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
    }
  } else if(scheme == MGGA_TAU) {
    switch(sel) {
      case DEN_S:
        zmat_mgga_vxc_onedft_kernel<false, DEN_S><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
      case DEN_Z:
        zmat_mgga_vxc_onedft_kernel<false, DEN_Z><<< blocks, threads, 0, stream >>>( ntasks, tasks_device );
        break;
    }
  } else {
    GAUXC_GENERIC_EXCEPTION("ONEDFT NYI for this scheme");
  }
}

void sz_to_ab(  size_t sz,
  const void* src_a,
  const void* src_b,
  void* dest_a,
  void* dest_b,
  device_queue queue ){
  
  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads(cuda::warp_size,cuda::max_warps_per_thread_block,1);
  dim3 blocks( util::div_ceil( sz, threads.x ), 1, 1 );

  sz_to_ab<<<blocks, threads, 0, stream>>>(
    sz,
    static_cast<const double*>(src_a),
    static_cast<const double*>(src_b),
    static_cast<double*>(dest_a),
    static_cast<double*>(dest_b)
  );
}

} // namespace GauXC::detail