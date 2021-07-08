#pragma once
#include "collocation_device_constants.hpp"
#include <cassert>
#include <type_traits>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __noinline__
#endif

namespace GauXC      {
namespace integrator {
namespace cuda       {

template <typename T, typename U>
struct collocation_cartesian_angular_impl;

template <typename T>
struct collocation_cartesian_angular_impl<T, std::integral_constant<int,0>> {

static GPGAUEVAL_INLINE __device__ void flat(
  uint32_t npts,
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  eval[npts * 0] = bf;

}

static GPGAUEVAL_INLINE __device__ void deriv1(
  uint32_t npts,
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

  eval_x[npts * 0] = bf_x;

  eval_y[npts * 0] = bf_y;

  eval_z[npts * 0] = bf_z;

}

};

template <typename T>
struct collocation_cartesian_angular_impl<T, std::integral_constant<int,1>> {

static GPGAUEVAL_INLINE __device__ void flat(
  uint32_t npts,
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  eval[npts * 0] = bf*x;
  eval[npts * 1] = bf*y;
  eval[npts * 2] = bf*z;

}

static GPGAUEVAL_INLINE __device__ void deriv1(
  uint32_t npts,
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

  eval_x[npts * 0] = bf + bf_x*x;
  eval_x[npts * 1] = bf_x*y;
  eval_x[npts * 2] = bf_x*z;

  eval_y[npts * 0] = bf_y*x;
  eval_y[npts * 1] = bf + bf_y*y;
  eval_y[npts * 2] = bf_y*z;

  eval_z[npts * 0] = bf_z*x;
  eval_z[npts * 1] = bf_z*y;
  eval_z[npts * 2] = bf + bf_z*z;

}

};

template <typename T>
struct collocation_cartesian_angular_impl<T, std::integral_constant<int,2>> {

static GPGAUEVAL_INLINE __device__ void flat(
  uint32_t npts,
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  eval[npts * 0] = bf*x*x;
  eval[npts * 1] = bf*x*y;
  eval[npts * 2] = bf*x*z;
  eval[npts * 3] = bf*y*y;
  eval[npts * 4] = bf*y*z;
  eval[npts * 5] = bf*z*z;

}

static GPGAUEVAL_INLINE __device__ void deriv1(
  uint32_t npts,
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

  eval_x[npts * 0] = x*(2*bf + bf_x*x);
  eval_x[npts * 1] = y*(bf + bf_x*x);
  eval_x[npts * 2] = z*(bf + bf_x*x);
  eval_x[npts * 3] = bf_x*y*y;
  eval_x[npts * 4] = bf_x*y*z;
  eval_x[npts * 5] = bf_x*z*z;

  eval_y[npts * 0] = bf_y*x*x;
  eval_y[npts * 1] = x*(bf + bf_y*y);
  eval_y[npts * 2] = bf_y*x*z;
  eval_y[npts * 3] = y*(2*bf + bf_y*y);
  eval_y[npts * 4] = z*(bf + bf_y*y);
  eval_y[npts * 5] = bf_y*z*z;

  eval_z[npts * 0] = bf_z*x*x;
  eval_z[npts * 1] = bf_z*x*y;
  eval_z[npts * 2] = x*(bf + bf_z*z);
  eval_z[npts * 3] = bf_z*y*y;
  eval_z[npts * 4] = y*(bf + bf_z*z);
  eval_z[npts * 5] = z*(2*bf + bf_z*z);

}

};

template <typename T>
struct collocation_cartesian_angular_impl<T, std::integral_constant<int,3>> {

static GPGAUEVAL_INLINE __device__ void flat(
  uint32_t npts,
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  eval[npts * 0] = bf*x*x*x;
  eval[npts * 1] = bf*x*x*y;
  eval[npts * 2] = bf*x*x*z;
  eval[npts * 3] = bf*x*y*y;
  eval[npts * 4] = bf*x*y*z;
  eval[npts * 5] = bf*x*z*z;
  eval[npts * 6] = bf*y*y*y;
  eval[npts * 7] = bf*y*y*z;
  eval[npts * 8] = bf*y*z*z;
  eval[npts * 9] = bf*z*z*z;

}

static GPGAUEVAL_INLINE __device__ void deriv1(
  uint32_t npts,
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

  eval_x[npts * 0] = x*x*(3*bf + bf_x*x);
  eval_x[npts * 1] = x*y*(2*bf + bf_x*x);
  eval_x[npts * 2] = x*z*(2*bf + bf_x*x);
  eval_x[npts * 3] = y*y*(bf + bf_x*x);
  eval_x[npts * 4] = y*z*(bf + bf_x*x);
  eval_x[npts * 5] = z*z*(bf + bf_x*x);
  eval_x[npts * 6] = bf_x*y*y*y;
  eval_x[npts * 7] = bf_x*y*y*z;
  eval_x[npts * 8] = bf_x*y*z*z;
  eval_x[npts * 9] = bf_x*z*z*z;

  eval_y[npts * 0] = bf_y*x*x*x;
  eval_y[npts * 1] = x*x*(bf + bf_y*y);
  eval_y[npts * 2] = bf_y*x*x*z;
  eval_y[npts * 3] = x*y*(2*bf + bf_y*y);
  eval_y[npts * 4] = x*z*(bf + bf_y*y);
  eval_y[npts * 5] = bf_y*x*z*z;
  eval_y[npts * 6] = y*y*(3*bf + bf_y*y);
  eval_y[npts * 7] = y*z*(2*bf + bf_y*y);
  eval_y[npts * 8] = z*z*(bf + bf_y*y);
  eval_y[npts * 9] = bf_y*z*z*z;

  eval_z[npts * 0] = bf_z*x*x*x;
  eval_z[npts * 1] = bf_z*x*x*y;
  eval_z[npts * 2] = x*x*(bf + bf_z*z);
  eval_z[npts * 3] = bf_z*x*y*y;
  eval_z[npts * 4] = x*y*(bf + bf_z*z);
  eval_z[npts * 5] = x*z*(2*bf + bf_z*z);
  eval_z[npts * 6] = bf_z*y*y*y;
  eval_z[npts * 7] = y*y*(bf + bf_z*z);
  eval_z[npts * 8] = y*z*(2*bf + bf_z*z);
  eval_z[npts * 9] = z*z*(3*bf + bf_z*z);

}

};


template <typename T, int L>
GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular(
  uint32_t npts,
  const T  bf,
  const T  x,
  const T  y,
  const T  z,
  T*       eval
) {

  collocation_cartesian_angular_impl<T, std::integral_constant<int, L>>::
    flat( npts, bf, x, y, z, eval );

} // collocation_cartesian_angular


template <typename T, int L>
GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_deriv1(
  uint32_t npts,
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

  collocation_cartesian_angular_impl<T, std::integral_constant<int, L>>::
    flat( npts, bf, x, y, z, eval );
  collocation_cartesian_angular_impl<T, std::integral_constant<int, L>>::
    deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

} // collocation_cartesian_angular_deriv1




template <typename T>
GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular(
  uint32_t         npts,
  const int32_t    l,
  const T          bf,
  const T          x,
  const T          y,
  const T          z,
  T*  eval
) {

      if( l == 0 ) {
  
        collocation_cartesian_angular<T,0>( npts, bf, x, y, z, eval );

      } else if( l == 1 ) {
  
        collocation_cartesian_angular<T,1>( npts, bf, x, y, z, eval );

      } else if( l == 2 ) {
  
        collocation_cartesian_angular<T,2>( npts, bf, x, y, z, eval );

      } else if( l == 3 ) {
  
        collocation_cartesian_angular<T,3>( npts, bf, x, y, z, eval );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_cartesian_angular


template <typename T>
GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_deriv1(
  uint32_t         npts,
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


      if( l == 0 ) {
  

        collocation_cartesian_angular_deriv1<T,0>( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
  

        collocation_cartesian_angular_deriv1<T,1>( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
  

        collocation_cartesian_angular_deriv1<T,2>( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
  

        collocation_cartesian_angular_deriv1<T,3>( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval, eval_x, eval_y, eval_z );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_cartesian_angular_deriv1


} // namespace cuda
} // namespace integrator
} // namespace GauXC

