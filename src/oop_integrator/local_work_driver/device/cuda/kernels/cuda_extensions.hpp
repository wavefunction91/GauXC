#pragma once
#include <cuda.h>
#include <cub/cub.cuh>

namespace GauXC {
namespace cuda  {

template <size_t warp_size_use, size_t warp_per_thread_block, typename T>
__device__ T warp_reduce_sum(T val) {

  using warp_reducer = cub::WarpReduce<T>;
  static __shared__ typename warp_reducer::TempStorage temp_storage[warp_per_thread_block];
  int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
  int warp_lane = tid / warp_size_use;
  val = warp_reducer( temp_storage[warp_lane] ).Sum( val );

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

  using warp_reducer = cub::WarpReduce<T>;
  static __shared__ typename warp_reducer::TempStorage temp_storage[warp_per_thread_block];
  int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
  int warp_lane = tid / warp_size_use;
  val = warp_reducer( temp_storage[warp_lane] ).Reduce( val, prod_operation() );

  return val;
}

}
}
