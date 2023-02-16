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
  auto* basis_eval_device = task.bf;
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

  //printf("[GPU] ITASK %d TID %d WL %d WID %d\n", batch_idx, threadIdx.x, warp_lane, warp_id);

  // Compute Max Bfn Sum
#if 0
  double max_bfn_sum = 0.0;
  for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
    double tmp = 0.0;
    for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
      tmp += std::abs( basis_eval_device[ ipt + ibf*npts ] );
    }
    max_bfn_sum = fmax( max_bfn_sum, std::sqrt(weights_device[ipt]) * tmp );
  }
  max_bfn_sum = cuda::warp_reduce_max<cuda::warp_size>(max_bfn_sum);
#else


  // First scale the basis functions by the weights
  for(int ipt = warp_lane; ipt < npts; ipt += cuda::warp_size) {
    const auto w = std::sqrt(weights_device[ipt]);
    for(int ibf = warp_id; ibf < nbf; ibf += nwarp) {
      const auto val = basis_eval_device[ ipt + ibf*npts ];
      basis_eval_device[ ipt + ibf*npts ] = w * std::abs(val);
    }
  }
  __syncthreads();




  __shared__ double bf_shared[32][32];
  __shared__ double bfn_sum_shared[32];
  bfn_sum_shared[warp_lane] = 0.0;
  __syncthreads();

  const int npts_chunks = GauXC::util::div_ceil(npts, cuda::warp_size);
  const int  nbf_chunks = GauXC::util::div_ceil(nbf,  nwarp);
  for(int ipts_chunk = 0; ipts_chunk < npts_chunks; ++ipts_chunk) {
    double tmp_bfn_sum = 0.0;
    const int ipt = ipts_chunk * cuda::warp_size + warp_lane;
  for(int  ibf_chunk = 0;  ibf_chunk <  nbf_chunks; ++ ibf_chunk) {
    const int ibf =  ibf_chunk * nwarp + warp_id;

    bf_shared[warp_id][warp_lane] = 0.0;

    // Load in a block of basis functions
    // Warp lane is the point index and warp ID is the bfn idx
    if(ipt < npts and ibf < nbf) 
      bf_shared[warp_id][warp_lane]  = basis_eval_device[ipt + ibf*npts];
    __syncthreads();

    // Do transpose
    // Warp lane is the bfn index and warp ID is the point idx
    auto tmp = bf_shared[warp_lane][warp_id];
    __syncthreads();
 
    // Do a sum reduce over basis functions for the chunk
    // Now every warp has the local bfn chunk sum in lane 0
    // corresponding to the point represented by the warp Id
    tmp_bfn_sum += cuda::warp_reduce_sum<cuda::warp_size>(tmp);
    
  }
    // At this point, every warp contains the total bfn sum
    // for the point corresponding to its warp id. Update the 
    // local value accordingly
    if(warp_lane == 0) {
      double val = bfn_sum_shared[warp_id];
      bfn_sum_shared[warp_id] = fmax( val, tmp_bfn_sum );
    }
    __syncthreads();

  }

  // Get global maximum
  double max_bfn_sum;
  if(warp_id == 0) {
    auto tmp = bfn_sum_shared[warp_lane];
    max_bfn_sum = cuda::warp_reduce_max<cuda::warp_size>(tmp);
  }

#endif
  if(threadIdx.x == 0) {
    task.max_bfn_sum = max_bfn_sum;
    printf("[GPU] ITASK = %d MAX_SUM = %.6e\n", batch_idx, max_bfn_sum);
    //printf("[GPU] ITASK = %d NBE = %lu NPTS = %lu \n", batch_idx, nbf, npts);
  }


}



void exx_ek_screening_bfn_stats( size_t        ntasks,
                                 XCDeviceTask* tasks_device,
                                 device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>() ;
  dim3 threads = 1024;//cuda::max_threads_per_thread_block;
  dim3 blocks  = ntasks;
  exx_ek_screening_bfn_stats_kernel<<<blocks, threads, 0, stream >>>(
    ntasks, tasks_device );

}

}
