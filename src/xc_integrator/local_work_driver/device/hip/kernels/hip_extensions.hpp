#pragma once
#include "hip/hip_runtime.h"
#include <hipcub/hipcub.hpp>
#include "device_specific/hip_device_constants.hpp"

namespace GauXC {
namespace hip   {

template <size_t warp_sz, typename T>
__device__ T warp_reduce_sum( T val ) { 

  using warp_reducer = hipcub::WarpReduce<double>;
  static __shared__ typename warp_reducer::TempStorage 
    temp_storage[hip::max_warps_per_thread_block];
  int tid = 
    threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;

  int warp_lane = tid / warp_size;

  return warp_reducer( temp_storage[warp_lane] ).Sum( val );

}

}
}
