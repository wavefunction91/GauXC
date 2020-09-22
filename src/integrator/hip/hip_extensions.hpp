#include "hip/hip_runtime.h"
#pragma once
#include <hip/hip_runtime.h>

namespace GauXC {
namespace hip  {

__inline__ __device__
double warpReduceSum(double val) {
  for(int i=16; i>=1; i/=2)
    val += __shfl_xor_sync(0xffffffff, val, i, 32);

  return val;
}

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

}
}
