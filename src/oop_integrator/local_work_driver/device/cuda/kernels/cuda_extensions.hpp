#pragma once
#include <cuda.h>
#include <cub/cub.cuh>

namespace GauXC {
namespace cuda  {

template <size_t warp_size_use, size_t warp_per_thread_block, typename T>
__device__ T warp_reduce_sum(T val) {

#if 0
  using warp_reducer = cub::WarpReduce<T>;
  static __shared__ typename warp_reducer::TempStorage temp_storage[warp_per_thread_block];
  int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
  int warp_lane = tid / warp_size_use;
  val = warp_reducer( temp_storage[warp_lane] ).Sum( val );
#else
  for(int i=(warp_size_use/2); i>=1; i/=2)
    val += __shfl_xor_sync(0xffffffff, val, i, warp_size_use);
#endif

  return val;
}

struct prod_operation {
  template <typename T>
  __device__  T operator()(const T& a, const T& b) const {
    return a*b;
  }
};

template <size_t warp_size_use, size_t warp_per_thread_block, typename T>
__device__ T warp_reduce_prod(T val) {

#if 0
  using warp_reducer = cub::WarpReduce<T>;
  static __shared__ typename warp_reducer::TempStorage temp_storage[warp_per_thread_block];
  int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
  int warp_lane = tid / warp_size_use;
  val = warp_reducer( temp_storage[warp_lane] ).Reduce( val, prod_operation() );
#else
  for(int i=(warp_size_use/2); i>=1; i/=2)
    val *= __shfl_xor_sync(0xffffffff, val, i, warp_size_use);
#endif

  return val;
}

}
}
