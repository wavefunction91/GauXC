#pragma once
#include "collocation_device_constants.hpp"
#include <cassert>
#include <type_traits>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace cuda       {

template <typename T, typename U>
struct collocation_$(name)_angular_impl;

$for( L in range(L_max + 1) )\
template <typename T>
struct collocation_$(name)_angular_impl<T, std::integral_constant<int,$(L)>> {

static GPGAUEVAL_INLINE __device__ void flat(
  const T bf,
  const T x,
  const T y,
  const T z,
  T*      eval
) {

  $(body[L])

}

static GPGAUEVAL_INLINE __device__ void deriv1(
  const T bf,
  const T bf_x,
  const T bf_y,
  const T bf_z,
  const T x,
  const T y,
  const T z,
  T* eval_x,
  T* eval_y,
  T* eval_z
) {

  $(body_d1[L])

}

};

$endfor\

template <typename T, int L>
GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular(
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  collocation_$(name)_angular_impl<T, std::integral_constant<int, L>>::
    flat( bf, x, y, z, eval );

} // collocation_$(name)_angular


template <typename T, int L>
GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular_deriv1(
  const T  bf,
  const T  bf_x,
  const T  bf_y,
  const T  bf_z,
  const T  x,
  const T  y,
  const T  z,
  T*       eval,
  T*       eval_x,
  T*       eval_y,
  T*       eval_z
) {

  collocation_$(name)_angular_impl<T, std::integral_constant<int, L>>::
    flat( bf, x, y, z, eval );
  collocation_$(name)_angular_impl<T, std::integral_constant<int, L>>::
    deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

} // collocation_$(name)_angular_deriv1




template <typename T>
GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular(
  const int32_t    l,
  const T          bf,
  const T          x,
  const T          y,
  const T          z,
  T*  eval
) {

$for( L in range(L_max + 1) )\
  $if( L == 0 )\
    if( l == $(L) ) {
  $else\
    } else if( l == $(L) ) {
  $endif
        collocation_$(name)_angular<T,$(L)>( bf, x, y, z, eval );

$endfor\
    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_$(name)_angular


template <typename T>
GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular_deriv1(
  const int32_t    l,
  const T          bf,
  const T          bf_x,
  const T          bf_y,
  const T          bf_z,
  const T          x,
  const T          y,
  const T          z,
  T*  eval,
  T*  eval_x,
  T*  eval_y,
  T*  eval_z
) {


$for( L in range(L_max + 1) )\
  $if( L == 0 )\
    if( l == $(L) ) {
  $else\
    } else if( l == $(L) ) {
  $endif

        collocation_$(name)_angular_deriv1<T,$(L)>( bf, bf_x, bf_y, bf_z, x, y, z, eval, eval_x, eval_y, eval_z );

$endfor\
    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_$(name)_angular_deriv1


} // namespace cuda
} // namespace integrator
} // namespace GauXC

