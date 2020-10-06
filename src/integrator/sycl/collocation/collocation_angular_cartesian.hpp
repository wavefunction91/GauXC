#pragma once
#include <CL/sycl.hpp>

#include "collocation_device_constants.hpp"
#include <gauxc/exceptions/gauxc_exception.hpp>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC      {
namespace integrator {
namespace sycl       {

GPGAUEVAL_INLINE void collocation_cartesian_angular_0(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_0_deriv1(
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

GPGAUEVAL_INLINE void collocation_cartesian_angular_1(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x;
  eval[1] = bf*y;
  eval[2] = bf*z;

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_1_deriv1(
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

  eval_x[0] = bf + bf_x*x;
  eval_x[1] = bf_x*y;
  eval_x[2] = bf_x*z;

  eval_y[0] = bf_y*x;
  eval_y[1] = bf + bf_y*y;
  eval_y[2] = bf_y*z;

  eval_z[0] = bf_z*x;
  eval_z[1] = bf_z*y;
  eval_z[2] = bf + bf_z*z;

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_2(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*x;
  eval[1] = bf*x*y;
  eval[2] = bf*x*z;
  eval[3] = bf*y*y;
  eval[4] = bf*y*z;
  eval[5] = bf*z*z;

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_2_deriv1(
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

  eval_x[0] = x*(2*bf + bf_x*x);
  eval_x[1] = y*(bf + bf_x*x);
  eval_x[2] = z*(bf + bf_x*x);
  eval_x[3] = bf_x*y*y;
  eval_x[4] = bf_x*y*z;
  eval_x[5] = bf_x*z*z;

  eval_y[0] = bf_y*x*x;
  eval_y[1] = x*(bf + bf_y*y);
  eval_y[2] = bf_y*x*z;
  eval_y[3] = y*(2*bf + bf_y*y);
  eval_y[4] = z*(bf + bf_y*y);
  eval_y[5] = bf_y*z*z;

  eval_z[0] = bf_z*x*x;
  eval_z[1] = bf_z*x*y;
  eval_z[2] = x*(bf + bf_z*z);
  eval_z[3] = bf_z*y*y;
  eval_z[4] = y*(bf + bf_z*z);
  eval_z[5] = z*(2*bf + bf_z*z);

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_3(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*x*x;
  eval[1] = bf*x*x*y;
  eval[2] = bf*x*x*z;
  eval[3] = bf*x*y*y;
  eval[4] = bf*x*y*z;
  eval[5] = bf*x*z*z;
  eval[6] = bf*y*y*y;
  eval[7] = bf*y*y*z;
  eval[8] = bf*y*z*z;
  eval[9] = bf*z*z*z;

}

GPGAUEVAL_INLINE void collocation_cartesian_angular_3_deriv1(
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

  eval_x[0] = x*x*(3*bf + bf_x*x);
  eval_x[1] = x*y*(2*bf + bf_x*x);
  eval_x[2] = x*z*(2*bf + bf_x*x);
  eval_x[3] = y*y*(bf + bf_x*x);
  eval_x[4] = y*z*(bf + bf_x*x);
  eval_x[5] = z*z*(bf + bf_x*x);
  eval_x[6] = bf_x*y*y*y;
  eval_x[7] = bf_x*y*y*z;
  eval_x[8] = bf_x*y*z*z;
  eval_x[9] = bf_x*z*z*z;

  eval_y[0] = bf_y*x*x*x;
  eval_y[1] = x*x*(bf + bf_y*y);
  eval_y[2] = bf_y*x*x*z;
  eval_y[3] = x*y*(2*bf + bf_y*y);
  eval_y[4] = x*z*(bf + bf_y*y);
  eval_y[5] = bf_y*x*z*z;
  eval_y[6] = y*y*(3*bf + bf_y*y);
  eval_y[7] = y*z*(2*bf + bf_y*y);
  eval_y[8] = z*z*(bf + bf_y*y);
  eval_y[9] = bf_y*z*z*z;

  eval_z[0] = bf_z*x*x*x;
  eval_z[1] = bf_z*x*x*y;
  eval_z[2] = x*x*(bf + bf_z*z);
  eval_z[3] = bf_z*x*y*y;
  eval_z[4] = x*y*(bf + bf_z*z);
  eval_z[5] = x*z*(2*bf + bf_z*z);
  eval_z[6] = bf_z*y*y*y;
  eval_z[7] = y*y*(bf + bf_z*z);
  eval_z[8] = y*z*(2*bf + bf_z*z);
  eval_z[9] = z*z*(3*bf + bf_z*z);

}


GPGAUEVAL_INLINE void collocation_cartesian_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        collocation_cartesian_angular_0( bf, x, y, z, eval );

      } else if( l == 1 ) {
        collocation_cartesian_angular_1( bf, x, y, z, eval );

      } else if( l == 2 ) {
        collocation_cartesian_angular_2( bf, x, y, z, eval );

      } else if( l == 3 ) {
        collocation_cartesian_angular_3( bf, x, y, z, eval );

    } else {
      //GAUXC_BOOL_CHECK( "L < L_MAX", false );
    }

} // collocation_cartesian_angular


GPGAUEVAL_INLINE void collocation_cartesian_angular_deriv1(
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
        collocation_cartesian_angular_0( bf, x, y, z, eval );
      collocation_cartesian_angular_0_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        collocation_cartesian_angular_1( bf, x, y, z, eval );
      collocation_cartesian_angular_1_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        collocation_cartesian_angular_2( bf, x, y, z, eval );
      collocation_cartesian_angular_2_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        collocation_cartesian_angular_3( bf, x, y, z, eval );
      collocation_cartesian_angular_3_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      //GAUXC_BOOL_CHECK( "L < L_MAX", false );
    }

} // collocation_cartesian_angular_deriv1



} // namespace sycl
} // namespace integrator
} // namespace GauXC
