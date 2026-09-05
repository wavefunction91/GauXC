#include "hip/hip_runtime.h"
#pragma once
#include <hip/hip_runtime.h>
#include <hipcub/hipcub.hpp>
#include "device/hip/hip_device_properties.hpp"

#define GAUXC_ENABLE_WARP_REDUCTIONS

namespace GauXC {
namespace hip  {

__inline__ __device__
double warpReduceSum(double val) {
 
#ifdef GAUXC_ENABLE_WARP_REDUCTIONS

  for(int i=(warp_size/2); i>=1; i/=2)
    val += __shfl_xor_sync(0xffffffff, val, i, warp_size);

#else

  using warp_reducer = hipcub::WarpReduce<double>;
  static __shared__ typename warp_reducer::TempStorage temp_storage[max_warps_per_thread_block];
  int tid = threadIdx.x + threadIdx.y * blockDim.x + threadIdx.z * blockDim.x * blockDim.y;
  int warp_lane = tid / warp_size;
  val = warp_reducer( temp_storage[warp_lane] ).Sum( val );

#endif

  return val;
}

__inline__ __device__
double warpReduceProd(double val) {
  for(int i=(warp_size/2); i>=1; i/=2)
    val *= __shfl_xor_sync(0xffffffff, val, i, warp_size);
  return val;
}

#if 0
__inline__ __device__
double blockReduceSum( double val ) {

  static __shared__ double shared[32];
  int lane = threadIdx.x % 32;
  int wid  = threadIdx.x / 32;

  val = warpReduceSum( val );

  if( lane == 0 ) shared[wid] = val;

  __syncthreads();

  val = (threadIdx.x < blockDim.x / 32) ? shared[lane] : 0;
  if( wid == 0 ) val = warpReduceSum( val );

  return val;

}

template <typename T, int warp_size = 32>
__inline__ __device__ T warp_prod_reduce( T val ) { 

  for( int i = warp_size / 2; i >= 1; i /= 2 )
    val *= __shfl_xor_sync( 0xffffffff, val, i, warp_size );

  return val;

}

template <typename T, int warp_size = 32 >
__inline__ __device__ T block_prod_reduce( T val ) {

  static __shared__ T shared[32];
  const int lane = threadIdx.x % 32;
  const int wid  = threadIdx.x / 32;

  val = warp_prod_reduce( val );

  if( lane == 0 ) shared[ wid ] = val;
  __syncthreads();

  val = ( threadIdx.x < blockDim.x / 32 ) ? shared[ lane ] : 0;
  if( wid == 0 ) val = warp_prod_reduce( val );

  return val;

}

__inline__ __device__ double atomicMul(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val *
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

}
}
