/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "device/common/exx_ek_screening.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"
#include "cuda_extensions.hpp"
#include "device_specific/cuda_device_constants.hpp"

namespace GauXC {

__global__ void exx_ek_screening_bfn_stats_kernel( size_t ntasks, 
                                                   XCDeviceTask* tasks_device ) {

  const int batch_idx = blockIdx.x;
  if( batch_idx >= ntasks ) return;
  
  auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  const auto* basis_eval_device = task.bf;
  const auto* weights_device    = task.weights;
  //double* bfn_max_device  = task.bfn_max;

  const int warp_lane = threadIdx.x % cuda::warp_size;
  const int warp_id   = threadIdx.x / cuda::warp_size;
  const int nwarp  = blockDim.x  / cuda::warp_size;

#if 0
  for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
    double tmp = 0.0;
    for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
      tmp = fmax( tmp,
        std::sqrt(weights_device[ipt]) *
        std::abs( basis_eval_device[ ipt + ibf*npts ] )
      );
    }
    bfn_max_device[ibf] = tmp;
  }
#endif

  // Compute Max Bfn Sum
  double max_bfn_sum = 0.0;
  for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
    double tmp = 0.0;
    for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
      tmp += std::abs( basis_eval_device[ ipt + ibf*npts ] );
    }
    max_bfn_sum = fmax( max_bfn_sum, std::sqrt(weights_device[ipt] * tmp) );
  }
  max_bfn_sum = cuda::warp_reduce_max<cuda::warp_size>(max_bfn_sum);
  if(threadIdx.x == 0) task.max_bfn_sum = max_bfn_sum;


}



void exx_ek_screening_bfn_stats( size_t        ntasks,
                                 XCDeviceTask* tasks_device,
                                 device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads = cuda::max_threads_per_thread_block;
  dim3 blocks  = ntasks;
  exx_ek_screening_bfn_stats_kernel<<<blocks, threads, 0, stream >>>(
    ntasks, tasks_device );

}

}
