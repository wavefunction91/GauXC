#pragma once
#include "collocation_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace cuda       {

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_0(
  const int npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0 * npts] = bf;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_0_deriv1(
  const int npts,
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

  eval_x[0 * npts] = bf_x;

  eval_y[0 * npts] = bf_y;

  eval_z[0 * npts] = bf_z;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_1(
  const int npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0 * npts] = bf*x;
  eval[1 * npts] = bf*y;
  eval[2 * npts] = bf*z;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_1_deriv1(
  const int npts,
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

  eval_x[0 * npts] = bf + bf_x*x;
  eval_x[1 * npts] = bf_x*y;
  eval_x[2 * npts] = bf_x*z;

  eval_y[0 * npts] = bf_y*x;
  eval_y[1 * npts] = bf + bf_y*y;
  eval_y[2 * npts] = bf_y*z;

  eval_z[0 * npts] = bf_z*x;
  eval_z[1 * npts] = bf_z*y;
  eval_z[2 * npts] = bf + bf_z*z;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_2(
  const int npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0 * npts] = bf*x*x;
  eval[1 * npts] = bf*x*y;
  eval[2 * npts] = bf*x*z;
  eval[3 * npts] = bf*y*y;
  eval[4 * npts] = bf*y*z;
  eval[5 * npts] = bf*z*z;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_2_deriv1(
  const int npts,
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

  eval_x[0 * npts] = x*(2*bf + bf_x*x);
  eval_x[1 * npts] = y*(bf + bf_x*x);
  eval_x[2 * npts] = z*(bf + bf_x*x);
  eval_x[3 * npts] = bf_x*y*y;
  eval_x[4 * npts] = bf_x*y*z;
  eval_x[5 * npts] = bf_x*z*z;

  eval_y[0 * npts] = bf_y*x*x;
  eval_y[1 * npts] = x*(bf + bf_y*y);
  eval_y[2 * npts] = bf_y*x*z;
  eval_y[3 * npts] = y*(2*bf + bf_y*y);
  eval_y[4 * npts] = z*(bf + bf_y*y);
  eval_y[5 * npts] = bf_y*z*z;

  eval_z[0 * npts] = bf_z*x*x;
  eval_z[1 * npts] = bf_z*x*y;
  eval_z[2 * npts] = x*(bf + bf_z*z);
  eval_z[3 * npts] = bf_z*y*y;
  eval_z[4 * npts] = y*(bf + bf_z*z);
  eval_z[5 * npts] = z*(2*bf + bf_z*z);

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_3(
  const int npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0 * npts] = bf*x*x*x;
  eval[1 * npts] = bf*x*x*y;
  eval[2 * npts] = bf*x*x*z;
  eval[3 * npts] = bf*x*y*y;
  eval[4 * npts] = bf*x*y*z;
  eval[5 * npts] = bf*x*z*z;
  eval[6 * npts] = bf*y*y*y;
  eval[7 * npts] = bf*y*y*z;
  eval[8 * npts] = bf*y*z*z;
  eval[9 * npts] = bf*z*z*z;

}

GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_3_deriv1(
  const int npts,
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

  eval_x[0 * npts] = x*x*(3*bf + bf_x*x);
  eval_x[1 * npts] = x*y*(2*bf + bf_x*x);
  eval_x[2 * npts] = x*z*(2*bf + bf_x*x);
  eval_x[3 * npts] = y*y*(bf + bf_x*x);
  eval_x[4 * npts] = y*z*(bf + bf_x*x);
  eval_x[5 * npts] = z*z*(bf + bf_x*x);
  eval_x[6 * npts] = bf_x*y*y*y;
  eval_x[7 * npts] = bf_x*y*y*z;
  eval_x[8 * npts] = bf_x*y*z*z;
  eval_x[9 * npts] = bf_x*z*z*z;

  eval_y[0 * npts] = bf_y*x*x*x;
  eval_y[1 * npts] = x*x*(bf + bf_y*y);
  eval_y[2 * npts] = bf_y*x*x*z;
  eval_y[3 * npts] = x*y*(2*bf + bf_y*y);
  eval_y[4 * npts] = x*z*(bf + bf_y*y);
  eval_y[5 * npts] = bf_y*x*z*z;
  eval_y[6 * npts] = y*y*(3*bf + bf_y*y);
  eval_y[7 * npts] = y*z*(2*bf + bf_y*y);
  eval_y[8 * npts] = z*z*(bf + bf_y*y);
  eval_y[9 * npts] = bf_y*z*z*z;

  eval_z[0 * npts] = bf_z*x*x*x;
  eval_z[1 * npts] = bf_z*x*x*y;
  eval_z[2 * npts] = x*x*(bf + bf_z*z);
  eval_z[3 * npts] = bf_z*x*y*y;
  eval_z[4 * npts] = x*y*(bf + bf_z*z);
  eval_z[5 * npts] = x*z*(2*bf + bf_z*z);
  eval_z[6 * npts] = bf_z*y*y*y;
  eval_z[7 * npts] = y*y*(bf + bf_z*z);
  eval_z[8 * npts] = y*z*(2*bf + bf_z*z);
  eval_z[9 * npts] = z*z*(3*bf + bf_z*z);

}


GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular(
  const int npts,
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        collocation_cartesian_angular_0( npts, bf, x, y, z, eval );

      } else if( l == 1 ) {
        collocation_cartesian_angular_1( npts, bf, x, y, z, eval );

      } else if( l == 2 ) {
        collocation_cartesian_angular_2( npts, bf, x, y, z, eval );

      } else if( l == 3 ) {
        collocation_cartesian_angular_3( npts, bf, x, y, z, eval );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_cartesian_angular


GPGAUEVAL_INLINE __device__ void collocation_cartesian_angular_deriv1(
  const int npts,
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


      if( l == 0 ) {
        collocation_cartesian_angular_0( npts, bf, x, y, z, eval );
      collocation_cartesian_angular_0_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        collocation_cartesian_angular_1( npts, bf, x, y, z, eval );
      collocation_cartesian_angular_1_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        collocation_cartesian_angular_2( npts, bf, x, y, z, eval );
      collocation_cartesian_angular_2_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        collocation_cartesian_angular_3( npts, bf, x, y, z, eval );
      collocation_cartesian_angular_3_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_cartesian_angular_deriv1



} // namespace cuda
} // namespace integrator
} // namespace GauXC
