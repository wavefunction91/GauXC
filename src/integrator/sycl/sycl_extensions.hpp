#pragma once
#include <SYCL/sycl.hpp>

namespace GauXC {
namespace sycl  {

__inline__
double warpReduceSum(double val, cl::sycl::nd_item<3> item_ct) {
  for(int i=16; i>=1; i/=2)
    val += item_ct.get_sub_group().shuffle_xor(val, i);

  return val;
}

__inline__
double blockReduceSum( double val , cl::sycl::nd_item<3> item_ct, double *shared) {

  int lane = item_ct.get_local_id(2) % 32;
  int wid = item_ct.get_local_id(2) / 32;

  val = warpReduceSum(val, item_ct);

  if( lane == 0 ) shared[wid] = val;

  group_barrier(item_ct.get_group());

  val = (item_ct.get_local_id(2) < item_ct.get_local_range().get(2) / 32)
            ? shared[lane]
            : 0;
    if (wid == 0) val = warpReduceSum(val, item_ct);

  return val;

}

template <typename T, int warp_size = 32>
__inline__ T warp_prod_reduce( cl::sycl::sub_group sub_g, T val ) {

  for( int i = warp_size / 2; i >= 1; i /= 2 )
      val *= sub_g.shuffle_xor( val, i );

  return val;
}

template <typename T, int warp_size = 32 >
__inline__ T block_prod_reduce( T val , cl::sycl::nd_item<3> item_ct, T *shared) {

  const int lane = item_ct.get_local_id(2) % 32;
  const int wid = item_ct.get_local_id(2) / 32;
  cl::sycl::sub_group sub_g = item_ct.get_sub_group();

  val = warp_prod_reduce( sub_g, val );

  if( lane == 0 ) shared[ wid ] = val;
  group_barrier(item_ct.get_group());

  val = (item_ct.get_local_id(2) < item_ct.get_local_range().get(2) / 32)
            ? shared[lane]
            : 0;
  if( wid == 0 ) val = warp_prod_reduce( sub_g, val );

  return val;

}

// __inline__ double atomicMul(double* address, double val)
// {
//     unsigned long long int* address_as_ull =
//                               (unsigned long long int*)address;
//     unsigned long long int old = *address_as_ull, assumed;
//
//     do {
//         assumed = old;
//     old = dpct::atomic_compare_exchange_strong(
//         address_as_ull, assumed,
//         (unsigned long long)(dpct::bit_cast<double, long long>(
//             val * dpct::bit_cast<long long, double>(assumed))));
//
//     // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
//     } while (assumed != old);
//
//   return dpct::bit_cast<long long, double>(old);
// }

}
}
