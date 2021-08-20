#pragma once
#include <cuda.h>

namespace GauXC {
namespace cuda  {

template <size_t warp_sz, typename T>
__device__ T warp_reduce_sum(T val) {

  for(int i=(warp_sz/2); i>=1; i/=2)
    val += __shfl_xor_sync(0xffffffff, val, i, warp_sz);

  return val;
}

template <size_t warp_sz, typename T>
__device__ T warp_reduce_prod(T val) {

  for(int i=(warp_sz/2); i>=1; i/=2)
    val *= __shfl_xor_sync(0xffffffff, val, i, warp_sz);

  return val;
}

}
}
