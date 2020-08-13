#pragma once
#include "collocation_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace cuda       {

$for( L in range(L_max + 1) )\
GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular_$(L)(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  $(body[L])

}

GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular_$(L)_deriv1(
  const double bf,
  const double bf_x,
  const double bf_y,
  const double bf_z,
  const double x,
  const double y,
  const double z,
  double* eval_x,
  double* eval_y,
  double* eval_z
) {

  $(body_d1[L])

}

$endfor\

GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

$for( L in range(L_max + 1) )\
  $if( L == 0 )\
    if( l == $(L) ) {
  $else\
    } else if( l == $(L) ) {
  $endif\
      collocation_$(name)_angular_$(L)( bf, x, y, z, eval );

$endfor\
    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_$(name)_angular


GPGAUEVAL_INLINE __device__ void collocation_$(name)_angular_deriv1(
  const int64_t l,
  const double  bf,
  const double  bf_x,
  const double  bf_y,
  const double  bf_z,
  const double  x,
  const double  y,
  const double  z,
  double*       eval,
  double*       eval_x,
  double*       eval_y,
  double*       eval_z
) {


$for( L in range(L_max + 1) )\
  $if( L == 0 )\
    if( l == $(L) ) {
  $else\
    } else if( l == $(L) ) {
  $endif\
      collocation_$(name)_angular_$(L)( bf, x, y, z, eval );
      collocation_$(name)_angular_$(L)_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

$endfor\
    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_$(name)_angular_deriv1


} // namespace cuda
} // namespace integrator
} // namespace GauXC

