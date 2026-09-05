#pragma once

#include "collocation_device_constants.hpp"
#include <gauxc/exceptions/gauxc_exception.hpp>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace sycl       {

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_0(
  const size_t npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_0_deriv1(
  const size_t npts,
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

  eval_x[0] = bf_x;

  eval_y[0] = bf_y;

  eval_z[0] = bf_z;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_1(
  const size_t npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[npts * 0] = bf*y;
  eval[npts * 1] = bf*z;
  eval[npts * 2] = bf*x;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_1_deriv1(
  const size_t npts,
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

  eval_x[npts * 0] = bf_x*y;
  eval_x[npts * 1] = bf_x*z;
  eval_x[npts * 2] = bf + bf_x*x;

  eval_y[npts * 0] = bf + bf_y*y;
  eval_y[npts * 1] = bf_y*z;
  eval_y[npts * 2] = bf_y*x;

  eval_z[npts * 0] = bf_z*y;
  eval_z[npts * 1] = bf + bf_z*z;
  eval_z[npts * 2] = bf_z*x;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_2(
  const size_t npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[npts * 0] = sqrt_3*bf*x*y;
  eval[npts * 1] = sqrt_3*bf*y*z;
  eval[npts * 2] = bf*(-x*x - y*y + 2*z*z)/2;
  eval[npts * 3] = sqrt_3*bf*x*z;
  eval[npts * 4] = sqrt_3*bf*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_2_deriv1(
  const size_t npts,
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

  eval_x[npts * 0] = sqrt_3*y*(bf + bf_x*x);
  eval_x[npts * 1] = sqrt_3*bf_x*y*z;
  eval_x[npts * 2] = -bf*x - bf_x*(x*x + y*y - 2*z*z)/2;
  eval_x[npts * 3] = sqrt_3*z*(bf + bf_x*x);
  eval_x[npts * 4] = sqrt_3*(bf*x + bf_x*(x*x - y*y)/2);

  eval_y[npts * 0] = sqrt_3*x*(bf + bf_y*y);
  eval_y[npts * 1] = sqrt_3*z*(bf + bf_y*y);
  eval_y[npts * 2] = -bf*y - bf_y*(x*x + y*y - 2*z*z)/2;
  eval_y[npts * 3] = sqrt_3*bf_y*x*z;
  eval_y[npts * 4] = sqrt_3*(-bf*y + bf_y*(x*x - y*y)/2);

  eval_z[npts * 0] = sqrt_3*bf_z*x*y;
  eval_z[npts * 1] = sqrt_3*y*(bf + bf_z*z);
  eval_z[npts * 2] = 2*bf*z - bf_z*(x*x + y*y - 2*z*z)/2;
  eval_z[npts * 3] = sqrt_3*x*(bf + bf_z*z);
  eval_z[npts * 4] = sqrt_3*bf_z*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_3(
  const size_t npts,
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[npts * 0] = sqrt_10*bf*y*(3*x*x - y*y)/4;
  eval[npts * 1] = sqrt_15*bf*x*y*z;
  eval[npts * 2] = sqrt_6*bf*y*(-x*x - y*y + 4*z*z)/4;
  eval[npts * 3] = bf*z*(-3*x*x - 3*y*y + 2*z*z)/2;
  eval[npts * 4] = sqrt_6*bf*x*(-x*x - y*y + 4*z*z)/4;
  eval[npts * 5] = sqrt_15*bf*z*(x*x - y*y)/2;
  eval[npts * 6] = sqrt_10*bf*x*(x*x - 3*y*y)/4;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_3_deriv1(
  const size_t npts,
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

  eval_x[npts * 0] = sqrt_10*y*(6*bf*x + bf_x*(3*x*x - y*y))/4;
  eval_x[npts * 1] = sqrt_15*y*z*(bf + bf_x*x);
  eval_x[npts * 2] = -sqrt_6*y*(2*bf*x + bf_x*(x*x + y*y - 4*z*z))/4;
  eval_x[npts * 3] = -z*(6*bf*x + bf_x*(3*x*x + 3*y*y - 2*z*z))/2;
  eval_x[npts * 4] = -sqrt_6*(bf*(3*x*x + y*y - 4*z*z) + bf_x*x*(x*x + y*y - 4*z*z))/4;
  eval_x[npts * 5] = sqrt_15*z*(2*bf*x + bf_x*(x*x - y*y))/2;
  eval_x[npts * 6] = sqrt_10*(3*bf*(x*x - y*y) + bf_x*x*(x*x - 3*y*y))/4;

  eval_y[npts * 0] = sqrt_10*(-3*bf*(-x*x + y*y) + bf_y*y*(3*x*x - y*y))/4;
  eval_y[npts * 1] = sqrt_15*x*z*(bf + bf_y*y);
  eval_y[npts * 2] = -sqrt_6*(bf*(x*x + 3*y*y - 4*z*z) + bf_y*y*(x*x + y*y - 4*z*z))/4;
  eval_y[npts * 3] = -z*(6*bf*y + bf_y*(3*x*x + 3*y*y - 2*z*z))/2;
  eval_y[npts * 4] = -sqrt_6*x*(2*bf*y + bf_y*(x*x + y*y - 4*z*z))/4;
  eval_y[npts * 5] = sqrt_15*z*(-2*bf*y + bf_y*(x*x - y*y))/2;
  eval_y[npts * 6] = sqrt_10*x*(-6*bf*y + bf_y*(x*x - 3*y*y))/4;

  eval_z[npts * 0] = sqrt_10*bf_z*y*(3*x*x - y*y)/4;
  eval_z[npts * 1] = sqrt_15*x*y*(bf + bf_z*z);
  eval_z[npts * 2] = sqrt_6*y*(8*bf*z - bf_z*(x*x + y*y - 4*z*z))/4;
  eval_z[npts * 3] = -3*bf*(x*x + y*y - 2*z*z)/2 - bf_z*z*(3*x*x + 3*y*y - 2*z*z)/2;
  eval_z[npts * 4] = sqrt_6*x*(8*bf*z - bf_z*(x*x + y*y - 4*z*z))/4;
  eval_z[npts * 5] = sqrt_15*(bf + bf_z*z)*(x*x - y*y)/2;
  eval_z[npts * 6] = sqrt_10*bf_z*x*(x*x - 3*y*y)/4;

}


GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular(
  const size_t npts,
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        collocation_spherical_unnorm_angular_0( npts, bf, x, y, z, eval );

      } else if( l == 1 ) {
        collocation_spherical_unnorm_angular_1( npts, bf, x, y, z, eval );

      } else if( l == 2 ) {
        collocation_spherical_unnorm_angular_2( npts, bf, x, y, z, eval );

      } else if( l == 3 ) {
        collocation_spherical_unnorm_angular_3( npts, bf, x, y, z, eval );

    } else {
      //GAUXC_BOOL_CHECK( "L < L_MAX", false );
    }

} // collocation_spherical_unnorm_angular


GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_deriv1(
  const size_t npts,
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
        collocation_spherical_unnorm_angular_0( npts, bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_0_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        collocation_spherical_unnorm_angular_1( npts, bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_1_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        collocation_spherical_unnorm_angular_2( npts, bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_2_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        collocation_spherical_unnorm_angular_3( npts, bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_3_deriv1( npts, bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      //GAUXC_BOOL_CHECK( "L < L_MAX", false );
    }

} // collocation_spherical_unnorm_angular_deriv1



} // namespace sycl
} // namespace integrator
} // namespace GauXC
