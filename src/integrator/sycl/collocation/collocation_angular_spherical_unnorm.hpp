#pragma once

#include "collocation_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace sycl       {

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_0(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_0_deriv1(
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
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*y;
  eval[1] = bf*z;
  eval[2] = bf*x;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_1_deriv1(
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

  eval_x[0] = bf_x*y;
  eval_x[1] = bf_x*z;
  eval_x[2] = bf + bf_x*x;

  eval_y[0] = bf + bf_y*y;
  eval_y[1] = bf_y*z;
  eval_y[2] = bf_y*x;

  eval_z[0] = bf_z*y;
  eval_z[1] = bf + bf_z*z;
  eval_z[2] = bf_z*x;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_2(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = sqrt_3*bf*x*y;
  eval[1] = sqrt_3*bf*y*z;
  eval[2] = bf*(-x*x - y*y + 2*z*z)/2;
  eval[3] = sqrt_3*bf*x*z;
  eval[4] = sqrt_3*bf*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_2_deriv1(
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

  eval_x[0] = sqrt_3*y*(bf + bf_x*x);
  eval_x[1] = sqrt_3*bf_x*y*z;
  eval_x[2] = -bf*x - bf_x*(x*x + y*y - 2*z*z)/2;
  eval_x[3] = sqrt_3*z*(bf + bf_x*x);
  eval_x[4] = sqrt_3*(bf*x + bf_x*(x*x - y*y)/2);

  eval_y[0] = sqrt_3*x*(bf + bf_y*y);
  eval_y[1] = sqrt_3*z*(bf + bf_y*y);
  eval_y[2] = -bf*y - bf_y*(x*x + y*y - 2*z*z)/2;
  eval_y[3] = sqrt_3*bf_y*x*z;
  eval_y[4] = sqrt_3*(-bf*y + bf_y*(x*x - y*y)/2);

  eval_z[0] = sqrt_3*bf_z*x*y;
  eval_z[1] = sqrt_3*y*(bf + bf_z*z);
  eval_z[2] = 2*bf*z - bf_z*(x*x + y*y - 2*z*z)/2;
  eval_z[3] = sqrt_3*x*(bf + bf_z*z);
  eval_z[4] = sqrt_3*bf_z*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_3(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = sqrt_10*bf*y*(3*x*x - y*y)/4;
  eval[1] = sqrt_15*bf*x*y*z;
  eval[2] = sqrt_6*bf*y*(-x*x - y*y + 4*z*z)/4;
  eval[3] = bf*z*(-3*x*x - 3*y*y + 2*z*z)/2;
  eval[4] = sqrt_6*bf*x*(-x*x - y*y + 4*z*z)/4;
  eval[5] = sqrt_15*bf*z*(x*x - y*y)/2;
  eval[6] = sqrt_10*bf*x*(x*x - 3*y*y)/4;

}

GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_3_deriv1(
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

  eval_x[0] = sqrt_10*y*(6*bf*x + bf_x*(3*x*x - y*y))/4;
  eval_x[1] = sqrt_15*y*z*(bf + bf_x*x);
  eval_x[2] = -sqrt_6*y*(2*bf*x + bf_x*(x*x + y*y - 4*z*z))/4;
  eval_x[3] = -z*(6*bf*x + bf_x*(3*x*x + 3*y*y - 2*z*z))/2;
  eval_x[4] = -sqrt_6*(bf*(3*x*x + y*y - 4*z*z) + bf_x*x*(x*x + y*y - 4*z*z))/4;
  eval_x[5] = sqrt_15*z*(2*bf*x + bf_x*(x*x - y*y))/2;
  eval_x[6] = sqrt_10*(3*bf*(x*x - y*y) + bf_x*x*(x*x - 3*y*y))/4;

  eval_y[0] = sqrt_10*(-3*bf*(-x*x + y*y) + bf_y*y*(3*x*x - y*y))/4;
  eval_y[1] = sqrt_15*x*z*(bf + bf_y*y);
  eval_y[2] = -sqrt_6*(bf*(x*x + 3*y*y - 4*z*z) + bf_y*y*(x*x + y*y - 4*z*z))/4;
  eval_y[3] = -z*(6*bf*y + bf_y*(3*x*x + 3*y*y - 2*z*z))/2;
  eval_y[4] = -sqrt_6*x*(2*bf*y + bf_y*(x*x + y*y - 4*z*z))/4;
  eval_y[5] = sqrt_15*z*(-2*bf*y + bf_y*(x*x - y*y))/2;
  eval_y[6] = sqrt_10*x*(-6*bf*y + bf_y*(x*x - 3*y*y))/4;

  eval_z[0] = sqrt_10*bf_z*y*(3*x*x - y*y)/4;
  eval_z[1] = sqrt_15*x*y*(bf + bf_z*z);
  eval_z[2] = sqrt_6*y*(8*bf*z - bf_z*(x*x + y*y - 4*z*z))/4;
  eval_z[3] = -3*bf*(x*x + y*y - 2*z*z)/2 - bf_z*z*(3*x*x + 3*y*y - 2*z*z)/2;
  eval_z[4] = sqrt_6*x*(8*bf*z - bf_z*(x*x + y*y - 4*z*z))/4;
  eval_z[5] = sqrt_15*(bf + bf_z*z)*(x*x - y*y)/2;
  eval_z[6] = sqrt_10*bf_z*x*(x*x - 3*y*y)/4;

}


GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        collocation_spherical_unnorm_angular_0( bf, x, y, z, eval );

      } else if( l == 1 ) {
        collocation_spherical_unnorm_angular_1( bf, x, y, z, eval );

      } else if( l == 2 ) {
        collocation_spherical_unnorm_angular_2( bf, x, y, z, eval );

      } else if( l == 3 ) {
        collocation_spherical_unnorm_angular_3( bf, x, y, z, eval );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_spherical_unnorm_angular


GPGAUEVAL_INLINE void collocation_spherical_unnorm_angular_deriv1(
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
        collocation_spherical_unnorm_angular_0( bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_0_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        collocation_spherical_unnorm_angular_1( bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_1_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        collocation_spherical_unnorm_angular_2( bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_2_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        collocation_spherical_unnorm_angular_3( bf, x, y, z, eval );
      collocation_spherical_unnorm_angular_3_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      assert( false && "L < L_MAX" );
    }

} // collocation_spherical_unnorm_angular_deriv1



} // namespace sycl
} // namespace integrator
} // namespace GauXC
