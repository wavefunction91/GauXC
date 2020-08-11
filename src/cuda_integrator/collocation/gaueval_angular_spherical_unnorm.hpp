#pragma once
#include "gaueval_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace GauXC {

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_0(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_0_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_1_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_2(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_2_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_3(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_3_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_4(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = sqrt_35*bf*x*y*(x*x - y*y)/2;
  eval[1] = sqrt_70*bf*y*z*(3*x*x - y*y)/4;
  eval[2] = sqrt_5*bf*x*y*(-x*x - y*y + 6*z*z)/2;
  eval[3] = sqrt_10*bf*y*z*(-3*x*x - 3*y*y + 4*z*z)/4;
  eval[4] = bf*(3*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z + 8*z*z*z*z)/8;
  eval[5] = sqrt_10*bf*x*z*(-3*x*x - 3*y*y + 4*z*z)/4;
  eval[6] = sqrt_5*bf*(-x*x*x*x + 6*x*x*z*z + y*y*y*y - 6*y*y*z*z)/4;
  eval[7] = sqrt_70*bf*x*z*(x*x - 3*y*y)/4;
  eval[8] = sqrt_35*bf*(x*x*x*x - 6*x*x*y*y + y*y*y*y)/8;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_4_deriv1(
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

  eval_x[0] = sqrt_35*y*(bf*(3*x*x - y*y) + bf_x*x*(x*x - y*y))/2;
  eval_x[1] = sqrt_70*y*z*(6*bf*x + bf_x*(3*x*x - y*y))/4;
  eval_x[2] = -sqrt_5*y*(bf*(3*x*x + y*y - 6*z*z) + bf_x*x*(x*x + y*y - 6*z*z))/2;
  eval_x[3] = -sqrt_10*y*z*(6*bf*x + bf_x*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_x[4] = 3*bf*x*(x*x + y*y - 4*z*z)/2 + bf_x*(3*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z + 8*z*z*z*z)/8;
  eval_x[5] = -sqrt_10*z*(bf*(9*x*x + 3*y*y - 4*z*z) + bf_x*x*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_x[6] = -sqrt_5*(bf*x*(x*x - 3*z*z) + bf_x*(x*x*x*x - 6*x*x*z*z - y*y*y*y + 6*y*y*z*z)/4);
  eval_x[7] = sqrt_70*z*(3*bf*(x*x - y*y) + bf_x*x*(x*x - 3*y*y))/4;
  eval_x[8] = sqrt_35*(4*bf*x*(x*x - 3*y*y) + bf_x*(x*x*x*x - 6*x*x*y*y + y*y*y*y))/8;

  eval_y[0] = sqrt_35*x*(-bf*(-x*x + 3*y*y) + bf_y*y*(x*x - y*y))/2;
  eval_y[1] = sqrt_70*z*(-3*bf*(-x*x + y*y) + bf_y*y*(3*x*x - y*y))/4;
  eval_y[2] = -sqrt_5*x*(bf*(x*x + 3*y*y - 6*z*z) + bf_y*y*(x*x + y*y - 6*z*z))/2;
  eval_y[3] = -sqrt_10*z*(bf*(3*x*x + 9*y*y - 4*z*z) + bf_y*y*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_y[4] = 3*bf*y*(x*x + y*y - 4*z*z)/2 + bf_y*(3*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z + 8*z*z*z*z)/8;
  eval_y[5] = -sqrt_10*x*z*(6*bf*y + bf_y*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_y[6] = sqrt_5*(bf*y*(y*y - 3*z*z) - bf_y*(x*x*x*x - 6*x*x*z*z - y*y*y*y + 6*y*y*z*z)/4);
  eval_y[7] = sqrt_70*x*z*(-6*bf*y + bf_y*(x*x - 3*y*y))/4;
  eval_y[8] = sqrt_35*(-4*bf*y*(3*x*x - y*y) + bf_y*(x*x*x*x - 6*x*x*y*y + y*y*y*y))/8;

  eval_z[0] = sqrt_35*bf_z*x*y*(x*x - y*y)/2;
  eval_z[1] = sqrt_70*y*(bf + bf_z*z)*(3*x*x - y*y)/4;
  eval_z[2] = sqrt_5*x*y*(12*bf*z - bf_z*(x*x + y*y - 6*z*z))/2;
  eval_z[3] = sqrt_10*y*(3*bf*(-x*x - y*y + 4*z*z) - bf_z*z*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_z[4] = -2*bf*z*(3*x*x + 3*y*y - 2*z*z) + bf_z*(3*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z + 8*z*z*z*z)/8;
  eval_z[5] = sqrt_10*x*(3*bf*(-x*x - y*y + 4*z*z) - bf_z*z*(3*x*x + 3*y*y - 4*z*z))/4;
  eval_z[6] = sqrt_5*(12*bf*z*(x*x - y*y) - bf_z*(x*x*x*x - 6*x*x*z*z - y*y*y*y + 6*y*y*z*z))/4;
  eval_z[7] = sqrt_70*x*(bf + bf_z*z)*(x*x - 3*y*y)/4;
  eval_z[8] = sqrt_35*bf_z*(x*x*x*x - 6*x*x*y*y + y*y*y*y)/8;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_5(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = 3*sqrt_14*bf*y*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y)/16;
  eval[1] = 3*sqrt_35*bf*x*y*z*(x*x - y*y)/2;
  eval[2] = sqrt_70*bf*y*(-3*x*x*x*x - 2*x*x*y*y + 24*x*x*z*z + y*y*y*y - 8*y*y*z*z)/16;
  eval[3] = sqrt_105*bf*x*y*z*(-x*x - y*y + 2*z*z)/2;
  eval[4] = sqrt_15*bf*y*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z)/8;
  eval[5] = bf*z*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z)/8;
  eval[6] = sqrt_15*bf*x*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z)/8;
  eval[7] = sqrt_105*bf*z*(-x*x*x*x + 2*x*x*z*z + y*y*y*y - 2*y*y*z*z)/4;
  eval[8] = sqrt_70*bf*x*(-x*x*x*x + 2*x*x*y*y + 8*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z)/16;
  eval[9] = 3*sqrt_35*bf*z*(x*x*x*x - 6*x*x*y*y + y*y*y*y)/8;
  eval[10] = 3*sqrt_14*bf*x*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y)/16;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_5_deriv1(
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

  eval_x[0] = 3*sqrt_14*y*(20*bf*x*(x*x - y*y) + bf_x*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y))/16;
  eval_x[1] = 3*sqrt_35*y*z*(bf*(3*x*x - y*y) + bf_x*x*(x*x - y*y))/2;
  eval_x[2] = -sqrt_70*y*(4*bf*x*(3*x*x + y*y - 12*z*z) + bf_x*(3*x*x*x*x + 2*x*x*y*y - 24*x*x*z*z - y*y*y*y + 8*y*y*z*z))/16;
  eval_x[3] = -sqrt_105*y*z*(bf*(3*x*x + y*y - 2*z*z) + bf_x*x*(x*x + y*y - 2*z*z))/2;
  eval_x[4] = sqrt_15*y*(4*bf*x*(x*x + y*y - 6*z*z) + bf_x*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_x[5] = z*(20*bf*x*(3*x*x + 3*y*y - 4*z*z) + bf_x*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z))/8;
  eval_x[6] = sqrt_15*(bf*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + 4*x*x*(x*x + y*y - 6*z*z) + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z) + bf_x*x*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_x[7] = -sqrt_105*z*(4*bf*x*(x*x - z*z) + bf_x*(x*x*x*x - 2*x*x*z*z - y*y*y*y + 2*y*y*z*z))/4;
  eval_x[8] = sqrt_70*(bf*(-x*x*x*x + 2*x*x*y*y + 8*x*x*z*z + 4*x*x*(-x*x + y*y + 4*z*z) + 3*y*y*y*y - 24*y*y*z*z) + bf_x*x*(-x*x*x*x + 2*x*x*y*y + 8*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z))/16;
  eval_x[9] = 3*sqrt_35*z*(4*bf*x*(x*x - 3*y*y) + bf_x*(x*x*x*x - 6*x*x*y*y + y*y*y*y))/8;
  eval_x[10] = 3*sqrt_14*(5*bf*x*x*x*x - 30*bf*x*x*y*y + 5*bf*y*y*y*y + bf_x*x*x*x*x*x - 10*bf_x*x*x*x*y*y + 5*bf_x*x*y*y*y*y)/16;

  eval_y[0] = 3*sqrt_14*(5*bf*x*x*x*x - 30*bf*x*x*y*y + 5*bf*y*y*y*y + 5*bf_y*x*x*x*x*y - 10*bf_y*x*x*y*y*y + bf_y*y*y*y*y*y)/16;
  eval_y[1] = 3*sqrt_35*x*z*(-bf*(-x*x + 3*y*y) + bf_y*y*(x*x - y*y))/2;
  eval_y[2] = sqrt_70*(-3*bf*x*x*x*x - 6*bf*x*x*y*y + 24*bf*x*x*z*z + 5*bf*y*y*y*y - 24*bf*y*y*z*z - 3*bf_y*x*x*x*x*y - 2*bf_y*x*x*y*y*y + 24*bf_y*x*x*y*z*z + bf_y*y*y*y*y*y - 8*bf_y*y*y*y*z*z)/16;
  eval_y[3] = -sqrt_105*x*z*(bf*(x*x + 3*y*y - 2*z*z) + bf_y*y*(x*x + y*y - 2*z*z))/2;
  eval_y[4] = sqrt_15*(bf*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 4*y*y*(x*x + y*y - 6*z*z) + 8*z*z*z*z) + bf_y*y*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_y[5] = z*(20*bf*y*(3*x*x + 3*y*y - 4*z*z) + bf_y*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z))/8;
  eval_y[6] = sqrt_15*x*(4*bf*y*(x*x + y*y - 6*z*z) + bf_y*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_y[7] = sqrt_105*z*(4*bf*y*(y*y - z*z) - bf_y*(x*x*x*x - 2*x*x*z*z - y*y*y*y + 2*y*y*z*z))/4;
  eval_y[8] = sqrt_70*x*(4*bf*y*(x*x + 3*y*y - 12*z*z) + bf_y*(-x*x*x*x + 2*x*x*y*y + 8*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z))/16;
  eval_y[9] = 3*sqrt_35*z*(-4*bf*y*(3*x*x - y*y) + bf_y*(x*x*x*x - 6*x*x*y*y + y*y*y*y))/8;
  eval_y[10] = 3*sqrt_14*x*(-20*bf*y*(x*x - y*y) + bf_y*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y))/16;

  eval_z[0] = 3*sqrt_14*bf_z*y*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y)/16;
  eval_z[1] = 3*sqrt_35*x*y*(bf + bf_z*z)*(x*x - y*y)/2;
  eval_z[2] = sqrt_70*y*(16*bf*z*(3*x*x - y*y) - bf_z*(3*x*x*x*x + 2*x*x*y*y - 24*x*x*z*z - y*y*y*y + 8*y*y*z*z))/16;
  eval_z[3] = sqrt_105*x*y*(bf*(-x*x - y*y + 6*z*z) - bf_z*z*(x*x + y*y - 2*z*z))/2;
  eval_z[4] = sqrt_15*y*(-8*bf*z*(3*x*x + 3*y*y - 4*z*z) + bf_z*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_z[5] = bf*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z + z*z*(-80*x*x - 80*y*y + 32*z*z))/8 + bf_z*z*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z)/8;
  eval_z[6] = sqrt_15*x*(-8*bf*z*(3*x*x + 3*y*y - 4*z*z) + bf_z*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z))/8;
  eval_z[7] = sqrt_105*(bf*(-x*x*x*x + 2*x*x*z*z + y*y*y*y - 2*y*y*z*z + 4*z*z*(x*x - y*y)) - bf_z*z*(x*x*x*x - 2*x*x*z*z - y*y*y*y + 2*y*y*z*z))/4;
  eval_z[8] = sqrt_70*x*(16*bf*z*(x*x - 3*y*y) + bf_z*(-x*x*x*x + 2*x*x*y*y + 8*x*x*z*z + 3*y*y*y*y - 24*y*y*z*z))/16;
  eval_z[9] = 3*sqrt_35*(bf + bf_z*z)*(x*x*x*x - 6*x*x*y*y + y*y*y*y)/8;
  eval_z[10] = 3*sqrt_14*bf_z*x*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y)/16;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_6(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = sqrt_462*bf*x*y*(3*x*x*x*x - 10*x*x*y*y + 3*y*y*y*y)/16;
  eval[1] = 3*sqrt_154*bf*y*z*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y)/16;
  eval[2] = 3*sqrt_7*bf*x*y*(-x*x*x*x + 10*x*x*z*z + y*y*y*y - 10*y*y*z*z)/4;
  eval[3] = sqrt_210*bf*y*z*(-9*x*x*x*x - 6*x*x*y*y + 24*x*x*z*z + 3*y*y*y*y - 8*y*y*z*z)/16;
  eval[4] = sqrt_210*bf*x*y*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + y*y*y*y - 16*y*y*z*z + 16*z*z*z*z)/16;
  eval[5] = sqrt_21*bf*y*z*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z)/8;
  eval[6] = bf*(-5*x*x*x*x*x*x - 15*x*x*x*x*y*y + 90*x*x*x*x*z*z - 15*x*x*y*y*y*y + 180*x*x*y*y*z*z - 120*x*x*z*z*z*z - 5*y*y*y*y*y*y + 90*y*y*y*y*z*z - 120*y*y*z*z*z*z + 16*z*z*z*z*z*z)/16;
  eval[7] = sqrt_21*bf*x*z*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z)/8;
  eval[8] = sqrt_210*bf*(x*x*x*x*x*x + x*x*x*x*y*y - 16*x*x*x*x*z*z - x*x*y*y*y*y + 16*x*x*z*z*z*z - y*y*y*y*y*y + 16*y*y*y*y*z*z - 16*y*y*z*z*z*z)/32;
  eval[9] = sqrt_210*bf*x*z*(-3*x*x*x*x + 6*x*x*y*y + 8*x*x*z*z + 9*y*y*y*y - 24*y*y*z*z)/16;
  eval[10] = 3*sqrt_7*bf*(-x*x*x*x*x*x + 5*x*x*x*x*y*y + 10*x*x*x*x*z*z + 5*x*x*y*y*y*y - 60*x*x*y*y*z*z - y*y*y*y*y*y + 10*y*y*y*y*z*z)/16;
  eval[11] = 3*sqrt_154*bf*x*z*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y)/16;
  eval[12] = sqrt_462*bf*(x*x*x*x*x*x - 15*x*x*x*x*y*y + 15*x*x*y*y*y*y - y*y*y*y*y*y)/32;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_6_deriv1(
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

  eval_x[0] = sqrt_462*y*(15*bf*x*x*x*x - 30*bf*x*x*y*y + 3*bf*y*y*y*y + 3*bf_x*x*x*x*x*x - 10*bf_x*x*x*x*y*y + 3*bf_x*x*y*y*y*y)/16;
  eval_x[1] = 3*sqrt_154*y*z*(20*bf*x*(x*x - y*y) + bf_x*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y))/16;
  eval_x[2] = 3*sqrt_7*y*(-5*bf*x*x*x*x + 30*bf*x*x*z*z + bf*y*y*y*y - 10*bf*y*y*z*z - bf_x*x*x*x*x*x + 10*bf_x*x*x*x*z*z + bf_x*x*y*y*y*y - 10*bf_x*x*y*y*z*z)/4;
  eval_x[3] = -sqrt_210*y*z*(12*bf*x*(3*x*x + y*y - 4*z*z) + bf_x*(9*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z - 3*y*y*y*y + 8*y*y*z*z))/16;
  eval_x[4] = sqrt_210*y*(bf*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + 4*x*x*(x*x + y*y - 8*z*z) + y*y*y*y - 16*y*y*z*z + 16*z*z*z*z) + bf_x*x*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + y*y*y*y - 16*y*y*z*z + 16*z*z*z*z))/16;
  eval_x[5] = sqrt_21*y*z*(20*bf*x*(x*x + y*y - 2*z*z) + bf_x*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_x[6] = -15*bf*x*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z)/8 - bf_x*(5*x*x*x*x*x*x + 15*x*x*x*x*y*y - 90*x*x*x*x*z*z + 15*x*x*y*y*y*y - 180*x*x*y*y*z*z + 120*x*x*z*z*z*z + 5*y*y*y*y*y*y - 90*y*y*y*y*z*z + 120*y*y*z*z*z*z - 16*z*z*z*z*z*z)/16;
  eval_x[7] = sqrt_21*z*(bf*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 20*x*x*(x*x + y*y - 2*z*z) + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z) + bf_x*x*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_x[8] = sqrt_210*(2*bf*x*(3*x*x*x*x + 2*x*x*y*y - 32*x*x*z*z - y*y*y*y + 16*z*z*z*z) + bf_x*(x*x*x*x*x*x + x*x*x*x*y*y - 16*x*x*x*x*z*z - x*x*y*y*y*y + 16*x*x*z*z*z*z - y*y*y*y*y*y + 16*y*y*y*y*z*z - 16*y*y*z*z*z*z))/32;
  eval_x[9] = sqrt_210*z*(-15*bf*x*x*x*x + 18*bf*x*x*y*y + 24*bf*x*x*z*z + 9*bf*y*y*y*y - 24*bf*y*y*z*z - 3*bf_x*x*x*x*x*x + 6*bf_x*x*x*x*y*y + 8*bf_x*x*x*x*z*z + 9*bf_x*x*y*y*y*y - 24*bf_x*x*y*y*z*z)/16;
  eval_x[10] = 3*sqrt_7*(2*bf*x*(-3*x*x*x*x + 10*x*x*y*y + 20*x*x*z*z + 5*y*y*y*y - 60*y*y*z*z) + bf_x*(-x*x*x*x*x*x + 5*x*x*x*x*y*y + 10*x*x*x*x*z*z + 5*x*x*y*y*y*y - 60*x*x*y*y*z*z - y*y*y*y*y*y + 10*y*y*y*y*z*z))/16;
  eval_x[11] = 3*sqrt_154*z*(5*bf*x*x*x*x - 30*bf*x*x*y*y + 5*bf*y*y*y*y + bf_x*x*x*x*x*x - 10*bf_x*x*x*x*y*y + 5*bf_x*x*y*y*y*y)/16;
  eval_x[12] = sqrt_462*(6*bf*x*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y) + bf_x*(x*x*x*x*x*x - 15*x*x*x*x*y*y + 15*x*x*y*y*y*y - y*y*y*y*y*y))/32;

  eval_y[0] = sqrt_462*x*(3*bf*x*x*x*x - 30*bf*x*x*y*y + 15*bf*y*y*y*y + 3*bf_y*x*x*x*x*y - 10*bf_y*x*x*y*y*y + 3*bf_y*y*y*y*y*y)/16;
  eval_y[1] = 3*sqrt_154*z*(5*bf*x*x*x*x - 30*bf*x*x*y*y + 5*bf*y*y*y*y + 5*bf_y*x*x*x*x*y - 10*bf_y*x*x*y*y*y + bf_y*y*y*y*y*y)/16;
  eval_y[2] = 3*sqrt_7*x*(bf*(-x*x*x*x + 10*x*x*z*z + y*y*y*y - 10*y*y*z*z + 4*y*y*(y*y - 5*z*z)) - bf_y*y*(x*x*x*x - 10*x*x*z*z - y*y*y*y + 10*y*y*z*z))/4;
  eval_y[3] = sqrt_210*z*(-9*bf*x*x*x*x - 18*bf*x*x*y*y + 24*bf*x*x*z*z + 15*bf*y*y*y*y - 24*bf*y*y*z*z - 9*bf_y*x*x*x*x*y - 6*bf_y*x*x*y*y*y + 24*bf_y*x*x*y*z*z + 3*bf_y*y*y*y*y*y - 8*bf_y*y*y*y*z*z)/16;
  eval_y[4] = sqrt_210*x*(bf*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + y*y*y*y - 16*y*y*z*z + 4*y*y*(x*x + y*y - 8*z*z) + 16*z*z*z*z) + bf_y*y*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + y*y*y*y - 16*y*y*z*z + 16*z*z*z*z))/16;
  eval_y[5] = sqrt_21*z*(bf*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 20*y*y*(x*x + y*y - 2*z*z) + 8*z*z*z*z) + bf_y*y*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_y[6] = -15*bf*y*(x*x*x*x + 2*x*x*y*y - 12*x*x*z*z + y*y*y*y - 12*y*y*z*z + 8*z*z*z*z)/8 - bf_y*(5*x*x*x*x*x*x + 15*x*x*x*x*y*y - 90*x*x*x*x*z*z + 15*x*x*y*y*y*y - 180*x*x*y*y*z*z + 120*x*x*z*z*z*z + 5*y*y*y*y*y*y - 90*y*y*y*y*z*z + 120*y*y*z*z*z*z - 16*z*z*z*z*z*z)/16;
  eval_y[7] = sqrt_21*x*z*(20*bf*y*(x*x + y*y - 2*z*z) + bf_y*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_y[8] = sqrt_210*(-2*bf*y*(-x*x*x*x + 2*x*x*y*y + 3*y*y*y*y - 32*y*y*z*z + 16*z*z*z*z) + bf_y*(x*x*x*x*x*x + x*x*x*x*y*y - 16*x*x*x*x*z*z - x*x*y*y*y*y + 16*x*x*z*z*z*z - y*y*y*y*y*y + 16*y*y*y*y*z*z - 16*y*y*z*z*z*z))/32;
  eval_y[9] = sqrt_210*x*z*(12*bf*y*(x*x + 3*y*y - 4*z*z) + bf_y*(-3*x*x*x*x + 6*x*x*y*y + 8*x*x*z*z + 9*y*y*y*y - 24*y*y*z*z))/16;
  eval_y[10] = 3*sqrt_7*(2*bf*y*(5*x*x*x*x + 10*x*x*y*y - 60*x*x*z*z - 3*y*y*y*y + 20*y*y*z*z) + bf_y*(-x*x*x*x*x*x + 5*x*x*x*x*y*y + 10*x*x*x*x*z*z + 5*x*x*y*y*y*y - 60*x*x*y*y*z*z - y*y*y*y*y*y + 10*y*y*y*y*z*z))/16;
  eval_y[11] = 3*sqrt_154*x*z*(-20*bf*y*(x*x - y*y) + bf_y*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y))/16;
  eval_y[12] = sqrt_462*(-6*bf*y*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y) + bf_y*(x*x*x*x*x*x - 15*x*x*x*x*y*y + 15*x*x*y*y*y*y - y*y*y*y*y*y))/32;

  eval_z[0] = sqrt_462*bf_z*x*y*(3*x*x*x*x - 10*x*x*y*y + 3*y*y*y*y)/16;
  eval_z[1] = 3*sqrt_154*y*(bf + bf_z*z)*(5*x*x*x*x - 10*x*x*y*y + y*y*y*y)/16;
  eval_z[2] = 3*sqrt_7*x*y*(20*bf*z*(x*x - y*y) - bf_z*(x*x*x*x - 10*x*x*z*z - y*y*y*y + 10*y*y*z*z))/4;
  eval_z[3] = sqrt_210*y*(bf*(-9*x*x*x*x - 6*x*x*y*y + 24*x*x*z*z + 3*y*y*y*y - 8*y*y*z*z + 16*z*z*(3*x*x - y*y)) - bf_z*z*(9*x*x*x*x + 6*x*x*y*y - 24*x*x*z*z - 3*y*y*y*y + 8*y*y*z*z))/16;
  eval_z[4] = sqrt_210*x*y*(-32*bf*z*(x*x + y*y - 2*z*z) + bf_z*(x*x*x*x + 2*x*x*y*y - 16*x*x*z*z + y*y*y*y - 16*y*y*z*z + 16*z*z*z*z))/16;
  eval_z[5] = sqrt_21*y*(-bf*(-5*x*x*x*x - 10*x*x*y*y + 20*x*x*z*z - 5*y*y*y*y + 20*y*y*z*z - 8*z*z*z*z + z*z*(40*x*x + 40*y*y - 32*z*z)) + bf_z*z*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_z[6] = 3*bf*z*(15*x*x*x*x + 30*x*x*y*y - 40*x*x*z*z + 15*y*y*y*y - 40*y*y*z*z + 8*z*z*z*z)/4 - bf_z*(5*x*x*x*x*x*x + 15*x*x*x*x*y*y - 90*x*x*x*x*z*z + 15*x*x*y*y*y*y - 180*x*x*y*y*z*z + 120*x*x*z*z*z*z + 5*y*y*y*y*y*y - 90*y*y*y*y*z*z + 120*y*y*z*z*z*z - 16*z*z*z*z*z*z)/16;
  eval_z[7] = sqrt_21*x*(-bf*(-5*x*x*x*x - 10*x*x*y*y + 20*x*x*z*z - 5*y*y*y*y + 20*y*y*z*z - 8*z*z*z*z + z*z*(40*x*x + 40*y*y - 32*z*z)) + bf_z*z*(5*x*x*x*x + 10*x*x*y*y - 20*x*x*z*z + 5*y*y*y*y - 20*y*y*z*z + 8*z*z*z*z))/8;
  eval_z[8] = sqrt_210*(-bf*z*(x*x*x*x - 2*x*x*z*z - y*y*y*y + 2*y*y*z*z) + bf_z*(x*x*x*x*x*x + x*x*x*x*y*y - 16*x*x*x*x*z*z - x*x*y*y*y*y + 16*x*x*z*z*z*z - y*y*y*y*y*y + 16*y*y*y*y*z*z - 16*y*y*z*z*z*z)/32);
  eval_z[9] = sqrt_210*x*(bf*(-3*x*x*x*x + 6*x*x*y*y + 8*x*x*z*z + 9*y*y*y*y - 24*y*y*z*z + 16*z*z*(x*x - 3*y*y)) + bf_z*z*(-3*x*x*x*x + 6*x*x*y*y + 8*x*x*z*z + 9*y*y*y*y - 24*y*y*z*z))/16;
  eval_z[10] = 3*sqrt_7*(20*bf*z*(x*x*x*x - 6*x*x*y*y + y*y*y*y) + bf_z*(-x*x*x*x*x*x + 5*x*x*x*x*y*y + 10*x*x*x*x*z*z + 5*x*x*y*y*y*y - 60*x*x*y*y*z*z - y*y*y*y*y*y + 10*y*y*y*y*z*z))/16;
  eval_z[11] = 3*sqrt_154*x*(bf + bf_z*z)*(x*x*x*x - 10*x*x*y*y + 5*y*y*y*y)/16;
  eval_z[12] = sqrt_462*bf_z*(x*x*x*x*x*x - 15*x*x*x*x*y*y + 15*x*x*y*y*y*y - y*y*y*y*y*y)/32;

}


GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        gaueval_spherical_unnorm_angular_0( bf, x, y, z, eval );

      } else if( l == 1 ) {
        gaueval_spherical_unnorm_angular_1( bf, x, y, z, eval );

      } else if( l == 2 ) {
        gaueval_spherical_unnorm_angular_2( bf, x, y, z, eval );

      } else if( l == 3 ) {
        gaueval_spherical_unnorm_angular_3( bf, x, y, z, eval );

      } else if( l == 4 ) {
        gaueval_spherical_unnorm_angular_4( bf, x, y, z, eval );

      } else if( l == 5 ) {
        gaueval_spherical_unnorm_angular_5( bf, x, y, z, eval );

      } else if( l == 6 ) {
        gaueval_spherical_unnorm_angular_6( bf, x, y, z, eval );

    } else {
      assert( false && "L < L_MAX" );
    }

} // gaueval_spherical_unnorm_angular


GPGAUEVAL_INLINE __device__ void gaueval_spherical_unnorm_angular_deriv1(
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
        gaueval_spherical_unnorm_angular_0( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_0_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        gaueval_spherical_unnorm_angular_1( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_1_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        gaueval_spherical_unnorm_angular_2( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_2_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        gaueval_spherical_unnorm_angular_3( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_3_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 4 ) {
        gaueval_spherical_unnorm_angular_4( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_4_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 5 ) {
        gaueval_spherical_unnorm_angular_5( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_5_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 6 ) {
        gaueval_spherical_unnorm_angular_6( bf, x, y, z, eval );
      gaueval_spherical_unnorm_angular_6_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      assert( false && "L < L_MAX" );
    }

} // gaueval_spherical_unnorm_angular_deriv1



} // namespace GauXC
