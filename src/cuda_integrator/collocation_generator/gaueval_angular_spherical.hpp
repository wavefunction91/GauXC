#pragma once
#include "gaueval_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace gpgaueval {

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_0(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_0_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_1_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_2(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*y;
  eval[1] = bf*y*z;
  eval[2] = bf*(-x*x - y*y + 2*z*z)/2;
  eval[3] = bf*x*z;
  eval[4] = sqrt_3*bf*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_2_deriv1(
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

  eval_x[0] = y*(bf + bf_x*x);
  eval_x[1] = bf_x*y*z;
  eval_x[2] = -bf*x - bf_x*(x*x + y*y - 2*z*z)/2;
  eval_x[3] = z*(bf + bf_x*x);
  eval_x[4] = sqrt_3*(bf*x + bf_x*(x*x - y*y)/2);

  eval_y[0] = x*(bf + bf_y*y);
  eval_y[1] = z*(bf + bf_y*y);
  eval_y[2] = -bf*y - bf_y*(x*x + y*y - 2*z*z)/2;
  eval_y[3] = bf_y*x*z;
  eval_y[4] = sqrt_3*(-bf*y + bf_y*(x*x - y*y)/2);

  eval_z[0] = bf_z*x*y;
  eval_z[1] = y*(bf + bf_z*z);
  eval_z[2] = 2*bf*z - bf_z*(x*x + y*y - 2*z*z)/2;
  eval_z[3] = x*(bf + bf_z*z);
  eval_z[4] = sqrt_3*bf_z*(x*x - y*y)/2;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_3(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*y*(3*sqrt_2*x*x - sqrt_10*y*y)/4;
  eval[1] = bf*x*y*z;
  eval[2] = bf*y*(-sqrt_30*x*x - 5*sqrt_6*y*y + 4*sqrt_30*z*z)/20;
  eval[3] = bf*z*(-3*sqrt_5*x*x - 3*sqrt_5*y*y + 10*z*z)/10;
  eval[4] = bf*x*(-5*sqrt_6*x*x - sqrt_30*y*y + 4*sqrt_30*z*z)/20;
  eval[5] = sqrt_3*bf*z*(x*x - y*y)/2;
  eval[6] = bf*x*(sqrt_10*x*x - 3*sqrt_2*y*y)/4;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_3_deriv1(
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

  eval_x[0] = y*(6*sqrt_2*bf*x + bf_x*(3*sqrt_2*x*x - sqrt_10*y*y))/4;
  eval_x[1] = y*z*(bf + bf_x*x);
  eval_x[2] = -y*(2*sqrt_30*bf*x + bf_x*(sqrt_30*x*x + 5*sqrt_6*y*y - 4*sqrt_30*z*z))/20;
  eval_x[3] = -z*(6*sqrt_5*bf*x + bf_x*(3*sqrt_5*x*x + 3*sqrt_5*y*y - 10*z*z))/10;
  eval_x[4] = -bf*(15*sqrt_6*x*x + sqrt_30*y*y - 4*sqrt_30*z*z)/20 - bf_x*x*(5*sqrt_6*x*x + sqrt_30*y*y - 4*sqrt_30*z*z)/20;
  eval_x[5] = sqrt_3*z*(2*bf*x + bf_x*(x*x - y*y))/2;
  eval_x[6] = 3*bf*(sqrt_10*x*x - sqrt_2*y*y)/4 + bf_x*x*(sqrt_10*x*x - 3*sqrt_2*y*y)/4;

  eval_y[0] = 3*bf*(sqrt_2*x*x - sqrt_10*y*y)/4 + bf_y*y*(3*sqrt_2*x*x - sqrt_10*y*y)/4;
  eval_y[1] = x*z*(bf + bf_y*y);
  eval_y[2] = -bf*(sqrt_30*x*x + 15*sqrt_6*y*y - 4*sqrt_30*z*z)/20 - bf_y*y*(sqrt_30*x*x + 5*sqrt_6*y*y - 4*sqrt_30*z*z)/20;
  eval_y[3] = -z*(6*sqrt_5*bf*y + bf_y*(3*sqrt_5*x*x + 3*sqrt_5*y*y - 10*z*z))/10;
  eval_y[4] = -x*(2*sqrt_30*bf*y + bf_y*(5*sqrt_6*x*x + sqrt_30*y*y - 4*sqrt_30*z*z))/20;
  eval_y[5] = sqrt_3*z*(-2*bf*y + bf_y*(x*x - y*y))/2;
  eval_y[6] = x*(-6*sqrt_2*bf*y + bf_y*(sqrt_10*x*x - 3*sqrt_2*y*y))/4;

  eval_z[0] = bf_z*y*(3*sqrt_2*x*x - sqrt_10*y*y)/4;
  eval_z[1] = x*y*(bf + bf_z*z);
  eval_z[2] = y*(8*sqrt_30*bf*z - bf_z*(sqrt_30*x*x + 5*sqrt_6*y*y - 4*sqrt_30*z*z))/20;
  eval_z[3] = -3*bf*(sqrt_5*x*x + sqrt_5*y*y - 10*z*z)/10 - bf_z*z*(3*sqrt_5*x*x + 3*sqrt_5*y*y - 10*z*z)/10;
  eval_z[4] = x*(8*sqrt_30*bf*z - bf_z*(5*sqrt_6*x*x + sqrt_30*y*y - 4*sqrt_30*z*z))/20;
  eval_z[5] = sqrt_3*(bf + bf_z*z)*(x*x - y*y)/2;
  eval_z[6] = bf_z*x*(sqrt_10*x*x - 3*sqrt_2*y*y)/4;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_4(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = sqrt_5*bf*x*y*(x*x - y*y)/2;
  eval[1] = bf*y*z*(3*sqrt_2*x*x - sqrt_10*y*y)/4;
  eval[2] = bf*x*y*(-sqrt_35*x*x - sqrt_35*y*y + 6*sqrt_7*z*z)/14;
  eval[3] = bf*y*z*(-3*sqrt_14*x*x - 3*sqrt_70*y*y + 4*sqrt_70*z*z)/28;
  eval[4] = bf*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 24*sqrt_105*x*x*z*z + 105*y*y*y*y - 24*sqrt_105*y*y*z*z + 280*z*z*z*z)/280;
  eval[5] = bf*x*z*(-3*sqrt_70*x*x - 3*sqrt_14*y*y + 4*sqrt_70*z*z)/28;
  eval[6] = bf*(-7*sqrt_5*x*x*x*x + 6*sqrt_21*x*x*z*z + 7*sqrt_5*y*y*y*y - 6*sqrt_21*y*y*z*z)/28;
  eval[7] = bf*x*z*(sqrt_10*x*x - 3*sqrt_2*y*y)/4;
  eval[8] = bf*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_4_deriv1(
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

  eval_x[0] = sqrt_5*y*(bf*(3*x*x - y*y) + bf_x*x*(x*x - y*y))/2;
  eval_x[1] = y*z*(6*sqrt_2*bf*x + bf_x*(3*sqrt_2*x*x - sqrt_10*y*y))/4;
  eval_x[2] = -y*(bf*(3*sqrt_35*x*x + sqrt_35*y*y - 6*sqrt_7*z*z) + bf_x*x*(sqrt_35*x*x + sqrt_35*y*y - 6*sqrt_7*z*z))/14;
  eval_x[3] = -y*z*(6*sqrt_14*bf*x + bf_x*(3*sqrt_14*x*x + 3*sqrt_70*y*y - 4*sqrt_70*z*z))/28;
  eval_x[4] = 3*bf*x*(35*x*x + sqrt_105*y*y - 4*sqrt_105*z*z)/70 + bf_x*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 24*sqrt_105*x*x*z*z + 105*y*y*y*y - 24*sqrt_105*y*y*z*z + 280*z*z*z*z)/280;
  eval_x[5] = -z*(bf*(9*sqrt_70*x*x + 3*sqrt_14*y*y - 4*sqrt_70*z*z) + bf_x*x*(3*sqrt_70*x*x + 3*sqrt_14*y*y - 4*sqrt_70*z*z))/28;
  eval_x[6] = -bf*x*(7*sqrt_5*x*x - 3*sqrt_21*z*z)/7 - bf_x*(7*sqrt_5*x*x*x*x - 6*sqrt_21*x*x*z*z - 7*sqrt_5*y*y*y*y + 6*sqrt_21*y*y*z*z)/28;
  eval_x[7] = z*(3*bf*(sqrt_10*x*x - sqrt_2*y*y) + bf_x*x*(sqrt_10*x*x - 3*sqrt_2*y*y))/4;
  eval_x[8] = bf*x*(sqrt_35*x*x - 3*sqrt_3*y*y)/2 + bf_x*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;

  eval_y[0] = sqrt_5*x*(-bf*(-x*x + 3*y*y) + bf_y*y*(x*x - y*y))/2;
  eval_y[1] = z*(-3*bf*(-sqrt_2*x*x + sqrt_10*y*y) + bf_y*y*(3*sqrt_2*x*x - sqrt_10*y*y))/4;
  eval_y[2] = -x*(bf*(sqrt_35*x*x + 3*sqrt_35*y*y - 6*sqrt_7*z*z) + bf_y*y*(sqrt_35*x*x + sqrt_35*y*y - 6*sqrt_7*z*z))/14;
  eval_y[3] = -z*(bf*(3*sqrt_14*x*x + 9*sqrt_70*y*y - 4*sqrt_70*z*z) + bf_y*y*(3*sqrt_14*x*x + 3*sqrt_70*y*y - 4*sqrt_70*z*z))/28;
  eval_y[4] = 3*bf*y*(sqrt_105*x*x + 35*y*y - 4*sqrt_105*z*z)/70 + bf_y*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 24*sqrt_105*x*x*z*z + 105*y*y*y*y - 24*sqrt_105*y*y*z*z + 280*z*z*z*z)/280;
  eval_y[5] = -x*z*(6*sqrt_14*bf*y + bf_y*(3*sqrt_70*x*x + 3*sqrt_14*y*y - 4*sqrt_70*z*z))/28;
  eval_y[6] = bf*y*(7*sqrt_5*y*y - 3*sqrt_21*z*z)/7 - bf_y*(7*sqrt_5*x*x*x*x - 6*sqrt_21*x*x*z*z - 7*sqrt_5*y*y*y*y + 6*sqrt_21*y*y*z*z)/28;
  eval_y[7] = x*z*(-6*sqrt_2*bf*y + bf_y*(sqrt_10*x*x - 3*sqrt_2*y*y))/4;
  eval_y[8] = -bf*y*(3*sqrt_3*x*x - sqrt_35*y*y)/2 + bf_y*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;

  eval_z[0] = sqrt_5*bf_z*x*y*(x*x - y*y)/2;
  eval_z[1] = y*(bf + bf_z*z)*(3*sqrt_2*x*x - sqrt_10*y*y)/4;
  eval_z[2] = x*y*(12*sqrt_7*bf*z - bf_z*(sqrt_35*x*x + sqrt_35*y*y - 6*sqrt_7*z*z))/14;
  eval_z[3] = y*(3*bf*(-sqrt_14*x*x - sqrt_70*y*y + 4*sqrt_70*z*z) - bf_z*z*(3*sqrt_14*x*x + 3*sqrt_70*y*y - 4*sqrt_70*z*z))/28;
  eval_z[4] = -2*bf*z*(3*sqrt_105*x*x + 3*sqrt_105*y*y - 70*z*z)/35 + bf_z*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 24*sqrt_105*x*x*z*z + 105*y*y*y*y - 24*sqrt_105*y*y*z*z + 280*z*z*z*z)/280;
  eval_z[5] = x*(3*bf*(-sqrt_70*x*x - sqrt_14*y*y + 4*sqrt_70*z*z) - bf_z*z*(3*sqrt_70*x*x + 3*sqrt_14*y*y - 4*sqrt_70*z*z))/28;
  eval_z[6] = 3*sqrt_21*bf*z*(x*x - y*y)/7 - bf_z*(7*sqrt_5*x*x*x*x - 6*sqrt_21*x*x*z*z - 7*sqrt_5*y*y*y*y + 6*sqrt_21*y*y*z*z)/28;
  eval_z[7] = x*(bf + bf_z*z)*(sqrt_10*x*x - 3*sqrt_2*y*y)/4;
  eval_z[8] = bf_z*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_5(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*y*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y)/16;
  eval[1] = sqrt_5*bf*x*y*z*(x*x - y*y)/2;
  eval[2] = bf*y*(-3*sqrt_70*x*x*x*x - 2*sqrt_30*x*x*y*y + 24*sqrt_6*x*x*z*z + 3*sqrt_70*y*y*y*y - 8*sqrt_30*y*y*z*z)/48;
  eval[3] = sqrt_15*bf*x*y*z*(-x*x - y*y + 2*z*z)/6;
  eval[4] = bf*y*(7*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_7*x*x*z*z + 21*sqrt_15*y*y*y*y - 36*sqrt_35*y*y*z*z + 56*sqrt_15*z*z*z*z)/168;
  eval[5] = bf*z*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 40*sqrt_21*x*x*z*z + 105*y*y*y*y - 40*sqrt_21*y*y*z*z + 168*z*z*z*z)/168;
  eval[6] = bf*x*(21*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_35*x*x*z*z + 7*sqrt_15*y*y*y*y - 36*sqrt_7*y*y*z*z + 56*sqrt_15*z*z*z*z)/168;
  eval[7] = bf*z*(-sqrt_105*x*x*x*x + 6*sqrt_5*x*x*z*z + sqrt_105*y*y*y*y - 6*sqrt_5*y*y*z*z)/12;
  eval[8] = bf*x*(-3*sqrt_70*x*x*x*x + 2*sqrt_30*x*x*y*y + 8*sqrt_30*x*x*z*z + 3*sqrt_70*y*y*y*y - 24*sqrt_6*y*y*z*z)/48;
  eval[9] = bf*z*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;
  eval[10] = bf*x*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y)/16;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_5_deriv1(
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

  eval_x[0] = y*(20*bf*x*(sqrt_14*x*x - sqrt_6*y*y) + bf_x*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y))/16;
  eval_x[1] = sqrt_5*y*z*(bf*(3*x*x - y*y) + bf_x*x*(x*x - y*y))/2;
  eval_x[2] = -y*(4*bf*x*(3*sqrt_70*x*x + sqrt_30*y*y - 12*sqrt_6*z*z) + bf_x*(3*sqrt_70*x*x*x*x + 2*sqrt_30*x*x*y*y - 24*sqrt_6*x*x*z*z - 3*sqrt_70*y*y*y*y + 8*sqrt_30*y*y*z*z))/48;
  eval_x[3] = -sqrt_15*y*z*(bf*(3*x*x + y*y - 2*z*z) + bf_x*x*(x*x + y*y - 2*z*z))/6;
  eval_x[4] = y*(4*bf*x*(7*sqrt_15*x*x + 3*sqrt_35*y*y - 18*sqrt_7*z*z) + bf_x*(7*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_7*x*x*z*z + 21*sqrt_15*y*y*y*y - 36*sqrt_35*y*y*z*z + 56*sqrt_15*z*z*z*z))/168;
  eval_x[5] = z*(4*bf*x*(105*x*x + 3*sqrt_105*y*y - 20*sqrt_21*z*z) + bf_x*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 40*sqrt_21*x*x*z*z + 105*y*y*y*y - 40*sqrt_21*y*y*z*z + 168*z*z*z*z))/168;
  eval_x[6] = 5*sqrt_15*bf*x*x*x*x/8 + 3*sqrt_35*bf*x*x*y*y/28 - 9*sqrt_35*bf*x*x*z*z/14 + sqrt_15*bf*y*y*y*y/24 - 3*sqrt_7*bf*y*y*z*z/14 + sqrt_15*bf*z*z*z*z/3 + sqrt_15*bf_x*x*x*x*x*x/8 + sqrt_35*bf_x*x*x*x*y*y/28 - 3*sqrt_35*bf_x*x*x*x*z*z/14 + sqrt_15*bf_x*x*y*y*y*y/24 - 3*sqrt_7*bf_x*x*y*y*z*z/14 + sqrt_15*bf_x*x*z*z*z*z/3;
  eval_x[7] = -z*(4*bf*x*(sqrt_105*x*x - 3*sqrt_5*z*z) + bf_x*(sqrt_105*x*x*x*x - 6*sqrt_5*x*x*z*z - sqrt_105*y*y*y*y + 6*sqrt_5*y*y*z*z))/12;
  eval_x[8] = -5*sqrt_70*bf*x*x*x*x/16 + sqrt_30*bf*x*x*y*y/8 + sqrt_30*bf*x*x*z*z/2 + sqrt_70*bf*y*y*y*y/16 - sqrt_6*bf*y*y*z*z/2 - sqrt_70*bf_x*x*x*x*x*x/16 + sqrt_30*bf_x*x*x*x*y*y/24 + sqrt_30*bf_x*x*x*x*z*z/6 + sqrt_70*bf_x*x*y*y*y*y/16 - sqrt_6*bf_x*x*y*y*z*z/2;
  eval_x[9] = z*(4*bf*x*(sqrt_35*x*x - 3*sqrt_3*y*y) + bf_x*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y))/8;
  eval_x[10] = 15*sqrt_14*bf*x*x*x*x/16 - 15*sqrt_6*bf*x*x*y*y/8 + 5*sqrt_14*bf*y*y*y*y/16 + 3*sqrt_14*bf_x*x*x*x*x*x/16 - 5*sqrt_6*bf_x*x*x*x*y*y/8 + 5*sqrt_14*bf_x*x*y*y*y*y/16;

  eval_y[0] = 5*sqrt_14*bf*x*x*x*x/16 - 15*sqrt_6*bf*x*x*y*y/8 + 15*sqrt_14*bf*y*y*y*y/16 + 5*sqrt_14*bf_y*x*x*x*x*y/16 - 5*sqrt_6*bf_y*x*x*y*y*y/8 + 3*sqrt_14*bf_y*y*y*y*y*y/16;
  eval_y[1] = sqrt_5*x*z*(-bf*(-x*x + 3*y*y) + bf_y*y*(x*x - y*y))/2;
  eval_y[2] = -sqrt_70*bf*x*x*x*x/16 - sqrt_30*bf*x*x*y*y/8 + sqrt_6*bf*x*x*z*z/2 + 5*sqrt_70*bf*y*y*y*y/16 - sqrt_30*bf*y*y*z*z/2 - sqrt_70*bf_y*x*x*x*x*y/16 - sqrt_30*bf_y*x*x*y*y*y/24 + sqrt_6*bf_y*x*x*y*z*z/2 + sqrt_70*bf_y*y*y*y*y*y/16 - sqrt_30*bf_y*y*y*y*z*z/6;
  eval_y[3] = -sqrt_15*x*z*(bf*(x*x + 3*y*y - 2*z*z) + bf_y*y*(x*x + y*y - 2*z*z))/6;
  eval_y[4] = sqrt_15*bf*x*x*x*x/24 + 3*sqrt_35*bf*x*x*y*y/28 - 3*sqrt_7*bf*x*x*z*z/14 + 5*sqrt_15*bf*y*y*y*y/8 - 9*sqrt_35*bf*y*y*z*z/14 + sqrt_15*bf*z*z*z*z/3 + sqrt_15*bf_y*x*x*x*x*y/24 + sqrt_35*bf_y*x*x*y*y*y/28 - 3*sqrt_7*bf_y*x*x*y*z*z/14 + sqrt_15*bf_y*y*y*y*y*y/8 - 3*sqrt_35*bf_y*y*y*y*z*z/14 + sqrt_15*bf_y*y*z*z*z*z/3;
  eval_y[5] = z*(4*bf*y*(3*sqrt_105*x*x + 105*y*y - 20*sqrt_21*z*z) + bf_y*(105*x*x*x*x + 6*sqrt_105*x*x*y*y - 40*sqrt_21*x*x*z*z + 105*y*y*y*y - 40*sqrt_21*y*y*z*z + 168*z*z*z*z))/168;
  eval_y[6] = x*(4*bf*y*(3*sqrt_35*x*x + 7*sqrt_15*y*y - 18*sqrt_7*z*z) + bf_y*(21*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_35*x*x*z*z + 7*sqrt_15*y*y*y*y - 36*sqrt_7*y*y*z*z + 56*sqrt_15*z*z*z*z))/168;
  eval_y[7] = z*(4*bf*y*(sqrt_105*y*y - 3*sqrt_5*z*z) - bf_y*(sqrt_105*x*x*x*x - 6*sqrt_5*x*x*z*z - sqrt_105*y*y*y*y + 6*sqrt_5*y*y*z*z))/12;
  eval_y[8] = x*(4*bf*y*(sqrt_30*x*x + 3*sqrt_70*y*y - 12*sqrt_6*z*z) + bf_y*(-3*sqrt_70*x*x*x*x + 2*sqrt_30*x*x*y*y + 8*sqrt_30*x*x*z*z + 3*sqrt_70*y*y*y*y - 24*sqrt_6*y*y*z*z))/48;
  eval_y[9] = z*(-4*bf*y*(3*sqrt_3*x*x - sqrt_35*y*y) + bf_y*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y))/8;
  eval_y[10] = x*(-20*bf*y*(sqrt_6*x*x - sqrt_14*y*y) + bf_y*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y))/16;

  eval_z[0] = bf_z*y*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y)/16;
  eval_z[1] = sqrt_5*x*y*(bf + bf_z*z)*(x*x - y*y)/2;
  eval_z[2] = y*(16*bf*z*(3*sqrt_6*x*x - sqrt_30*y*y) - bf_z*(3*sqrt_70*x*x*x*x + 2*sqrt_30*x*x*y*y - 24*sqrt_6*x*x*z*z - 3*sqrt_70*y*y*y*y + 8*sqrt_30*y*y*z*z))/48;
  eval_z[3] = sqrt_15*x*y*(bf*(-x*x - y*y + 6*z*z) - bf_z*z*(x*x + y*y - 2*z*z))/6;
  eval_z[4] = y*(-8*bf*z*(9*sqrt_7*x*x + 9*sqrt_35*y*y - 28*sqrt_15*z*z) + bf_z*(7*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_7*x*x*z*z + 21*sqrt_15*y*y*y*y - 36*sqrt_35*y*y*z*z + 56*sqrt_15*z*z*z*z))/168;
  eval_z[5] = 5*bf*x*x*x*x/8 + sqrt_105*bf*x*x*y*y/28 - 5*sqrt_21*bf*x*x*z*z/7 + 5*bf*y*y*y*y/8 - 5*sqrt_21*bf*y*y*z*z/7 + 5*bf*z*z*z*z + 5*bf_z*x*x*x*x*z/8 + sqrt_105*bf_z*x*x*y*y*z/28 - 5*sqrt_21*bf_z*x*x*z*z*z/21 + 5*bf_z*y*y*y*y*z/8 - 5*sqrt_21*bf_z*y*y*z*z*z/21 + bf_z*z*z*z*z*z;
  eval_z[6] = x*(-8*bf*z*(9*sqrt_35*x*x + 9*sqrt_7*y*y - 28*sqrt_15*z*z) + bf_z*(21*sqrt_15*x*x*x*x + 6*sqrt_35*x*x*y*y - 36*sqrt_35*x*x*z*z + 7*sqrt_15*y*y*y*y - 36*sqrt_7*y*y*z*z + 56*sqrt_15*z*z*z*z))/168;
  eval_z[7] = -sqrt_105*bf*x*x*x*x/12 + 3*sqrt_5*bf*x*x*z*z/2 + sqrt_105*bf*y*y*y*y/12 - 3*sqrt_5*bf*y*y*z*z/2 - sqrt_105*bf_z*x*x*x*x*z/12 + sqrt_5*bf_z*x*x*z*z*z/2 + sqrt_105*bf_z*y*y*y*y*z/12 - sqrt_5*bf_z*y*y*z*z*z/2;
  eval_z[8] = x*(16*bf*z*(sqrt_30*x*x - 3*sqrt_6*y*y) + bf_z*(-3*sqrt_70*x*x*x*x + 2*sqrt_30*x*x*y*y + 8*sqrt_30*x*x*z*z + 3*sqrt_70*y*y*y*y - 24*sqrt_6*y*y*z*z))/48;
  eval_z[9] = (bf + bf_z*z)*(sqrt_35*x*x*x*x - 6*sqrt_3*x*x*y*y + sqrt_35*y*y*y*y)/8;
  eval_z[10] = bf_z*x*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y)/16;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_6(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*y*(3*sqrt_42*x*x*x*x - 10*sqrt_10*x*x*y*y + 3*sqrt_42*y*y*y*y)/16;
  eval[1] = bf*y*z*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y)/16;
  eval[2] = bf*x*y*(-3*sqrt_77*x*x*x*x + 10*sqrt_33*x*x*z*z + 3*sqrt_77*y*y*y*y - 10*sqrt_33*y*y*z*z)/44;
  eval[3] = bf*y*z*(-3*sqrt_2310*x*x*x*x - 6*sqrt_110*x*x*y*y + 24*sqrt_110*x*x*z*z + 3*sqrt_2310*y*y*y*y - 40*sqrt_22*y*y*z*z)/176;
  eval[4] = bf*x*y*(3*sqrt_2310*x*x*x*x + 30*sqrt_22*x*x*y*y - 48*sqrt_110*x*x*z*z + 3*sqrt_2310*y*y*y*y - 48*sqrt_110*y*y*z*z + 16*sqrt_2310*z*z*z*z)/528;
  eval[5] = bf*y*z*(5*sqrt_231*x*x*x*x + 30*sqrt_11*x*x*y*y - 60*sqrt_11*x*x*z*z + 15*sqrt_231*y*y*y*y - 60*sqrt_55*y*y*z*z + 24*sqrt_231*z*z*z*z)/264;
  eval[6] = bf*(-385*x*x*x*x*x*x - 35*sqrt_33*x*x*x*x*y*y + 210*sqrt_33*x*x*x*x*z*z - 35*sqrt_33*x*x*y*y*y*y + 36*sqrt_385*x*x*y*y*z*z - 280*sqrt_33*x*x*z*z*z*z - 385*y*y*y*y*y*y + 210*sqrt_33*y*y*y*y*z*z - 280*sqrt_33*y*y*z*z*z*z + 1232*z*z*z*z*z*z)/1232;
  eval[7] = bf*x*z*(15*sqrt_231*x*x*x*x + 30*sqrt_11*x*x*y*y - 60*sqrt_55*x*x*z*z + 5*sqrt_231*y*y*y*y - 60*sqrt_11*y*y*z*z + 24*sqrt_231*z*z*z*z)/264;
  eval[8] = bf*(11*sqrt_210*x*x*x*x*x*x + sqrt_770*x*x*x*x*y*y - 16*sqrt_770*x*x*x*x*z*z - sqrt_770*x*x*y*y*y*y + 16*sqrt_770*x*x*z*z*z*z - 11*sqrt_210*y*y*y*y*y*y + 16*sqrt_770*y*y*y*y*z*z - 16*sqrt_770*y*y*z*z*z*z)/352;
  eval[9] = bf*x*z*(-3*sqrt_2310*x*x*x*x + 6*sqrt_110*x*x*y*y + 40*sqrt_22*x*x*z*z + 3*sqrt_2310*y*y*y*y - 24*sqrt_110*y*y*z*z)/176;
  eval[10] = bf*(-33*sqrt_7*x*x*x*x*x*x + 5*sqrt_231*x*x*x*x*y*y + 10*sqrt_231*x*x*x*x*z*z + 5*sqrt_231*x*x*y*y*y*y - 36*sqrt_55*x*x*y*y*z*z - 33*sqrt_7*y*y*y*y*y*y + 10*sqrt_231*y*y*y*y*z*z)/176;
  eval[11] = bf*x*z*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y)/16;
  eval[12] = bf*(sqrt_462*x*x*x*x*x*x - 15*sqrt_14*x*x*x*x*y*y + 15*sqrt_14*x*x*y*y*y*y - sqrt_462*y*y*y*y*y*y)/32;

}

GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_6_deriv1(
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

  eval_x[0] = y*(15*sqrt_42*bf*x*x*x*x - 30*sqrt_10*bf*x*x*y*y + 3*sqrt_42*bf*y*y*y*y + 3*sqrt_42*bf_x*x*x*x*x*x - 10*sqrt_10*bf_x*x*x*x*y*y + 3*sqrt_42*bf_x*x*y*y*y*y)/16;
  eval_x[1] = y*z*(20*bf*x*(sqrt_14*x*x - sqrt_6*y*y) + bf_x*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y))/16;
  eval_x[2] = y*(-15*sqrt_77*bf*x*x*x*x + 30*sqrt_33*bf*x*x*z*z + 3*sqrt_77*bf*y*y*y*y - 10*sqrt_33*bf*y*y*z*z - 3*sqrt_77*bf_x*x*x*x*x*x + 10*sqrt_33*bf_x*x*x*x*z*z + 3*sqrt_77*bf_x*x*y*y*y*y - 10*sqrt_33*bf_x*x*y*y*z*z)/44;
  eval_x[3] = -y*z*(12*bf*x*(sqrt_2310*x*x + sqrt_110*y*y - 4*sqrt_110*z*z) + bf_x*(3*sqrt_2310*x*x*x*x + 6*sqrt_110*x*x*y*y - 24*sqrt_110*x*x*z*z - 3*sqrt_2310*y*y*y*y + 40*sqrt_22*y*y*z*z))/176;
  eval_x[4] = y*(15*sqrt_2310*bf*x*x*x*x + 90*sqrt_22*bf*x*x*y*y - 144*sqrt_110*bf*x*x*z*z + 3*sqrt_2310*bf*y*y*y*y - 48*sqrt_110*bf*y*y*z*z + 16*sqrt_2310*bf*z*z*z*z + 3*sqrt_2310*bf_x*x*x*x*x*x + 30*sqrt_22*bf_x*x*x*x*y*y - 48*sqrt_110*bf_x*x*x*x*z*z + 3*sqrt_2310*bf_x*x*y*y*y*y - 48*sqrt_110*bf_x*x*y*y*z*z + 16*sqrt_2310*bf_x*x*z*z*z*z)/528;
  eval_x[5] = y*z*(20*bf*x*(sqrt_231*x*x + 3*sqrt_11*y*y - 6*sqrt_11*z*z) + bf_x*(5*sqrt_231*x*x*x*x + 30*sqrt_11*x*x*y*y - 60*sqrt_11*x*x*z*z + 15*sqrt_231*y*y*y*y - 60*sqrt_55*y*y*z*z + 24*sqrt_231*z*z*z*z))/264;
  eval_x[6] = -bf*x*(1155*x*x*x*x + 70*sqrt_33*x*x*y*y - 420*sqrt_33*x*x*z*z + 35*sqrt_33*y*y*y*y - 36*sqrt_385*y*y*z*z + 280*sqrt_33*z*z*z*z)/616 - bf_x*(385*x*x*x*x*x*x + 35*sqrt_33*x*x*x*x*y*y - 210*sqrt_33*x*x*x*x*z*z + 35*sqrt_33*x*x*y*y*y*y - 36*sqrt_385*x*x*y*y*z*z + 280*sqrt_33*x*x*z*z*z*z + 385*y*y*y*y*y*y - 210*sqrt_33*y*y*y*y*z*z + 280*sqrt_33*y*y*z*z*z*z - 1232*z*z*z*z*z*z)/1232;
  eval_x[7] = z*(75*sqrt_231*bf*x*x*x*x + 90*sqrt_11*bf*x*x*y*y - 180*sqrt_55*bf*x*x*z*z + 5*sqrt_231*bf*y*y*y*y - 60*sqrt_11*bf*y*y*z*z + 24*sqrt_231*bf*z*z*z*z + 15*sqrt_231*bf_x*x*x*x*x*x + 30*sqrt_11*bf_x*x*x*x*y*y - 60*sqrt_55*bf_x*x*x*x*z*z + 5*sqrt_231*bf_x*x*y*y*y*y - 60*sqrt_11*bf_x*x*y*y*z*z + 24*sqrt_231*bf_x*x*z*z*z*z)/264;
  eval_x[8] = bf*x*(33*sqrt_210*x*x*x*x + 2*sqrt_770*x*x*y*y - 32*sqrt_770*x*x*z*z - sqrt_770*y*y*y*y + 16*sqrt_770*z*z*z*z)/176 + bf_x*(11*sqrt_210*x*x*x*x*x*x + sqrt_770*x*x*x*x*y*y - 16*sqrt_770*x*x*x*x*z*z - sqrt_770*x*x*y*y*y*y + 16*sqrt_770*x*x*z*z*z*z - 11*sqrt_210*y*y*y*y*y*y + 16*sqrt_770*y*y*y*y*z*z - 16*sqrt_770*y*y*z*z*z*z)/352;
  eval_x[9] = z*(-15*sqrt_2310*bf*x*x*x*x + 18*sqrt_110*bf*x*x*y*y + 120*sqrt_22*bf*x*x*z*z + 3*sqrt_2310*bf*y*y*y*y - 24*sqrt_110*bf*y*y*z*z - 3*sqrt_2310*bf_x*x*x*x*x*x + 6*sqrt_110*bf_x*x*x*x*y*y + 40*sqrt_22*bf_x*x*x*x*z*z + 3*sqrt_2310*bf_x*x*y*y*y*y - 24*sqrt_110*bf_x*x*y*y*z*z)/176;
  eval_x[10] = bf*x*(-99*sqrt_7*x*x*x*x + 10*sqrt_231*x*x*y*y + 20*sqrt_231*x*x*z*z + 5*sqrt_231*y*y*y*y - 36*sqrt_55*y*y*z*z)/88 + bf_x*(-33*sqrt_7*x*x*x*x*x*x + 5*sqrt_231*x*x*x*x*y*y + 10*sqrt_231*x*x*x*x*z*z + 5*sqrt_231*x*x*y*y*y*y - 36*sqrt_55*x*x*y*y*z*z - 33*sqrt_7*y*y*y*y*y*y + 10*sqrt_231*y*y*y*y*z*z)/176;
  eval_x[11] = z*(15*sqrt_14*bf*x*x*x*x - 30*sqrt_6*bf*x*x*y*y + 5*sqrt_14*bf*y*y*y*y + 3*sqrt_14*bf_x*x*x*x*x*x - 10*sqrt_6*bf_x*x*x*x*y*y + 5*sqrt_14*bf_x*x*y*y*y*y)/16;
  eval_x[12] = 3*bf*x*(sqrt_462*x*x*x*x - 10*sqrt_14*x*x*y*y + 5*sqrt_14*y*y*y*y)/16 + bf_x*(sqrt_462*x*x*x*x*x*x - 15*sqrt_14*x*x*x*x*y*y + 15*sqrt_14*x*x*y*y*y*y - sqrt_462*y*y*y*y*y*y)/32;

  eval_y[0] = x*(3*sqrt_42*bf*x*x*x*x - 30*sqrt_10*bf*x*x*y*y + 15*sqrt_42*bf*y*y*y*y + 3*sqrt_42*bf_y*x*x*x*x*y - 10*sqrt_10*bf_y*x*x*y*y*y + 3*sqrt_42*bf_y*y*y*y*y*y)/16;
  eval_y[1] = z*(5*sqrt_14*bf*x*x*x*x - 30*sqrt_6*bf*x*x*y*y + 15*sqrt_14*bf*y*y*y*y + 5*sqrt_14*bf_y*x*x*x*x*y - 10*sqrt_6*bf_y*x*x*y*y*y + 3*sqrt_14*bf_y*y*y*y*y*y)/16;
  eval_y[2] = x*(-3*sqrt_77*bf*x*x*x*x + 10*sqrt_33*bf*x*x*z*z + 15*sqrt_77*bf*y*y*y*y - 30*sqrt_33*bf*y*y*z*z - 3*sqrt_77*bf_y*x*x*x*x*y + 10*sqrt_33*bf_y*x*x*y*z*z + 3*sqrt_77*bf_y*y*y*y*y*y - 10*sqrt_33*bf_y*y*y*y*z*z)/44;
  eval_y[3] = z*(-3*sqrt_2310*bf*x*x*x*x - 18*sqrt_110*bf*x*x*y*y + 24*sqrt_110*bf*x*x*z*z + 15*sqrt_2310*bf*y*y*y*y - 120*sqrt_22*bf*y*y*z*z - 3*sqrt_2310*bf_y*x*x*x*x*y - 6*sqrt_110*bf_y*x*x*y*y*y + 24*sqrt_110*bf_y*x*x*y*z*z + 3*sqrt_2310*bf_y*y*y*y*y*y - 40*sqrt_22*bf_y*y*y*y*z*z)/176;
  eval_y[4] = x*(3*sqrt_2310*bf*x*x*x*x + 90*sqrt_22*bf*x*x*y*y - 48*sqrt_110*bf*x*x*z*z + 15*sqrt_2310*bf*y*y*y*y - 144*sqrt_110*bf*y*y*z*z + 16*sqrt_2310*bf*z*z*z*z + 3*sqrt_2310*bf_y*x*x*x*x*y + 30*sqrt_22*bf_y*x*x*y*y*y - 48*sqrt_110*bf_y*x*x*y*z*z + 3*sqrt_2310*bf_y*y*y*y*y*y - 48*sqrt_110*bf_y*y*y*y*z*z + 16*sqrt_2310*bf_y*y*z*z*z*z)/528;
  eval_y[5] = z*(5*sqrt_231*bf*x*x*x*x + 90*sqrt_11*bf*x*x*y*y - 60*sqrt_11*bf*x*x*z*z + 75*sqrt_231*bf*y*y*y*y - 180*sqrt_55*bf*y*y*z*z + 24*sqrt_231*bf*z*z*z*z + 5*sqrt_231*bf_y*x*x*x*x*y + 30*sqrt_11*bf_y*x*x*y*y*y - 60*sqrt_11*bf_y*x*x*y*z*z + 15*sqrt_231*bf_y*y*y*y*y*y - 60*sqrt_55*bf_y*y*y*y*z*z + 24*sqrt_231*bf_y*y*z*z*z*z)/264;
  eval_y[6] = -bf*y*(35*sqrt_33*x*x*x*x + 70*sqrt_33*x*x*y*y - 36*sqrt_385*x*x*z*z + 1155*y*y*y*y - 420*sqrt_33*y*y*z*z + 280*sqrt_33*z*z*z*z)/616 - bf_y*(385*x*x*x*x*x*x + 35*sqrt_33*x*x*x*x*y*y - 210*sqrt_33*x*x*x*x*z*z + 35*sqrt_33*x*x*y*y*y*y - 36*sqrt_385*x*x*y*y*z*z + 280*sqrt_33*x*x*z*z*z*z + 385*y*y*y*y*y*y - 210*sqrt_33*y*y*y*y*z*z + 280*sqrt_33*y*y*z*z*z*z - 1232*z*z*z*z*z*z)/1232;
  eval_y[7] = x*z*(20*bf*y*(3*sqrt_11*x*x + sqrt_231*y*y - 6*sqrt_11*z*z) + bf_y*(15*sqrt_231*x*x*x*x + 30*sqrt_11*x*x*y*y - 60*sqrt_55*x*x*z*z + 5*sqrt_231*y*y*y*y - 60*sqrt_11*y*y*z*z + 24*sqrt_231*z*z*z*z))/264;
  eval_y[8] = -bf*y*(-sqrt_770*x*x*x*x + 2*sqrt_770*x*x*y*y + 33*sqrt_210*y*y*y*y - 32*sqrt_770*y*y*z*z + 16*sqrt_770*z*z*z*z)/176 + bf_y*(11*sqrt_210*x*x*x*x*x*x + sqrt_770*x*x*x*x*y*y - 16*sqrt_770*x*x*x*x*z*z - sqrt_770*x*x*y*y*y*y + 16*sqrt_770*x*x*z*z*z*z - 11*sqrt_210*y*y*y*y*y*y + 16*sqrt_770*y*y*y*y*z*z - 16*sqrt_770*y*y*z*z*z*z)/352;
  eval_y[9] = x*z*(12*bf*y*(sqrt_110*x*x + sqrt_2310*y*y - 4*sqrt_110*z*z) + bf_y*(-3*sqrt_2310*x*x*x*x + 6*sqrt_110*x*x*y*y + 40*sqrt_22*x*x*z*z + 3*sqrt_2310*y*y*y*y - 24*sqrt_110*y*y*z*z))/176;
  eval_y[10] = bf*y*(5*sqrt_231*x*x*x*x + 10*sqrt_231*x*x*y*y - 36*sqrt_55*x*x*z*z - 99*sqrt_7*y*y*y*y + 20*sqrt_231*y*y*z*z)/88 + bf_y*(-33*sqrt_7*x*x*x*x*x*x + 5*sqrt_231*x*x*x*x*y*y + 10*sqrt_231*x*x*x*x*z*z + 5*sqrt_231*x*x*y*y*y*y - 36*sqrt_55*x*x*y*y*z*z - 33*sqrt_7*y*y*y*y*y*y + 10*sqrt_231*y*y*y*y*z*z)/176;
  eval_y[11] = x*z*(-20*bf*y*(sqrt_6*x*x - sqrt_14*y*y) + bf_y*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y))/16;
  eval_y[12] = -3*bf*y*(5*sqrt_14*x*x*x*x - 10*sqrt_14*x*x*y*y + sqrt_462*y*y*y*y)/16 + bf_y*(sqrt_462*x*x*x*x*x*x - 15*sqrt_14*x*x*x*x*y*y + 15*sqrt_14*x*x*y*y*y*y - sqrt_462*y*y*y*y*y*y)/32;

  eval_z[0] = bf_z*x*y*(3*sqrt_42*x*x*x*x - 10*sqrt_10*x*x*y*y + 3*sqrt_42*y*y*y*y)/16;
  eval_z[1] = y*(bf + bf_z*z)*(5*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 3*sqrt_14*y*y*y*y)/16;
  eval_z[2] = x*y*(20*sqrt_33*bf*z*(x*x - y*y) - bf_z*(3*sqrt_77*x*x*x*x - 10*sqrt_33*x*x*z*z - 3*sqrt_77*y*y*y*y + 10*sqrt_33*y*y*z*z))/44;
  eval_z[3] = y*(-3*sqrt_2310*bf*x*x*x*x - 6*sqrt_110*bf*x*x*y*y + 72*sqrt_110*bf*x*x*z*z + 3*sqrt_2310*bf*y*y*y*y - 120*sqrt_22*bf*y*y*z*z - 3*sqrt_2310*bf_z*x*x*x*x*z - 6*sqrt_110*bf_z*x*x*y*y*z + 24*sqrt_110*bf_z*x*x*z*z*z + 3*sqrt_2310*bf_z*y*y*y*y*z - 40*sqrt_22*bf_z*y*y*z*z*z)/176;
  eval_z[4] = x*y*(-32*bf*z*(3*sqrt_110*x*x + 3*sqrt_110*y*y - 2*sqrt_2310*z*z) + bf_z*(3*sqrt_2310*x*x*x*x + 30*sqrt_22*x*x*y*y - 48*sqrt_110*x*x*z*z + 3*sqrt_2310*y*y*y*y - 48*sqrt_110*y*y*z*z + 16*sqrt_2310*z*z*z*z))/528;
  eval_z[5] = y*(5*sqrt_231*bf*x*x*x*x + 30*sqrt_11*bf*x*x*y*y - 180*sqrt_11*bf*x*x*z*z + 15*sqrt_231*bf*y*y*y*y - 180*sqrt_55*bf*y*y*z*z + 120*sqrt_231*bf*z*z*z*z + 5*sqrt_231*bf_z*x*x*x*x*z + 30*sqrt_11*bf_z*x*x*y*y*z - 60*sqrt_11*bf_z*x*x*z*z*z + 15*sqrt_231*bf_z*y*y*y*y*z - 60*sqrt_55*bf_z*y*y*z*z*z + 24*sqrt_231*bf_z*z*z*z*z*z)/264;
  eval_z[6] = bf*z*(105*sqrt_33*x*x*x*x + 18*sqrt_385*x*x*y*y - 280*sqrt_33*x*x*z*z + 105*sqrt_33*y*y*y*y - 280*sqrt_33*y*y*z*z + 1848*z*z*z*z)/308 - bf_z*(385*x*x*x*x*x*x + 35*sqrt_33*x*x*x*x*y*y - 210*sqrt_33*x*x*x*x*z*z + 35*sqrt_33*x*x*y*y*y*y - 36*sqrt_385*x*x*y*y*z*z + 280*sqrt_33*x*x*z*z*z*z + 385*y*y*y*y*y*y - 210*sqrt_33*y*y*y*y*z*z + 280*sqrt_33*y*y*z*z*z*z - 1232*z*z*z*z*z*z)/1232;
  eval_z[7] = x*(15*sqrt_231*bf*x*x*x*x + 30*sqrt_11*bf*x*x*y*y - 180*sqrt_55*bf*x*x*z*z + 5*sqrt_231*bf*y*y*y*y - 180*sqrt_11*bf*y*y*z*z + 120*sqrt_231*bf*z*z*z*z + 15*sqrt_231*bf_z*x*x*x*x*z + 30*sqrt_11*bf_z*x*x*y*y*z - 60*sqrt_55*bf_z*x*x*z*z*z + 5*sqrt_231*bf_z*y*y*y*y*z - 60*sqrt_11*bf_z*y*y*z*z*z + 24*sqrt_231*bf_z*z*z*z*z*z)/264;
  eval_z[8] = -sqrt_770*bf*z*(x*x*x*x - 2*x*x*z*z - y*y*y*y + 2*y*y*z*z)/11 + bf_z*(11*sqrt_210*x*x*x*x*x*x + sqrt_770*x*x*x*x*y*y - 16*sqrt_770*x*x*x*x*z*z - sqrt_770*x*x*y*y*y*y + 16*sqrt_770*x*x*z*z*z*z - 11*sqrt_210*y*y*y*y*y*y + 16*sqrt_770*y*y*y*y*z*z - 16*sqrt_770*y*y*z*z*z*z)/352;
  eval_z[9] = x*(-3*sqrt_2310*bf*x*x*x*x + 6*sqrt_110*bf*x*x*y*y + 120*sqrt_22*bf*x*x*z*z + 3*sqrt_2310*bf*y*y*y*y - 72*sqrt_110*bf*y*y*z*z - 3*sqrt_2310*bf_z*x*x*x*x*z + 6*sqrt_110*bf_z*x*x*y*y*z + 40*sqrt_22*bf_z*x*x*z*z*z + 3*sqrt_2310*bf_z*y*y*y*y*z - 24*sqrt_110*bf_z*y*y*z*z*z)/176;
  eval_z[10] = bf*z*(5*sqrt_231*x*x*x*x - 18*sqrt_55*x*x*y*y + 5*sqrt_231*y*y*y*y)/44 + bf_z*(-33*sqrt_7*x*x*x*x*x*x + 5*sqrt_231*x*x*x*x*y*y + 10*sqrt_231*x*x*x*x*z*z + 5*sqrt_231*x*x*y*y*y*y - 36*sqrt_55*x*x*y*y*z*z - 33*sqrt_7*y*y*y*y*y*y + 10*sqrt_231*y*y*y*y*z*z)/176;
  eval_z[11] = x*(bf + bf_z*z)*(3*sqrt_14*x*x*x*x - 10*sqrt_6*x*x*y*y + 5*sqrt_14*y*y*y*y)/16;
  eval_z[12] = bf_z*(sqrt_462*x*x*x*x*x*x - 15*sqrt_14*x*x*x*x*y*y + 15*sqrt_14*x*x*y*y*y*y - sqrt_462*y*y*y*y*y*y)/32;

}


GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

//switch( l ) {

  //case 0:
      if( l == 0 ) {
        gaueval_spherical_angular_0( bf, x, y, z, eval );
      //break;

  //case 1:
      } else if( l == 1 ) {
        gaueval_spherical_angular_1( bf, x, y, z, eval );
      //break;

  //case 2:
      } else if( l == 2 ) {
        gaueval_spherical_angular_2( bf, x, y, z, eval );
      //break;

  //case 3:
      } else if( l == 3 ) {
        gaueval_spherical_angular_3( bf, x, y, z, eval );
      //break;

  //case 4:
      } else if( l == 4 ) {
        gaueval_spherical_angular_4( bf, x, y, z, eval );
      //break;

  //case 5:
      } else if( l == 5 ) {
        gaueval_spherical_angular_5( bf, x, y, z, eval );
      //break;

  //case 6:
      } else if( l == 6 ) {
        gaueval_spherical_angular_6( bf, x, y, z, eval );
      //break;

    //default: 
    } else {
      assert( false && "L < L_MAX" );
      //break;
    }

//} // switch(l)

} // gaueval_spherical_angular


GPGAUEVAL_INLINE __device__ void gaueval_spherical_angular_deriv1(
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

//switch( l ) {

//  case 0:
      if( l == 0 ) {
        gaueval_spherical_angular_0( bf, x, y, z, eval );
      gaueval_spherical_angular_0_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 1:
      } else if( l == 1 ) {
        gaueval_spherical_angular_1( bf, x, y, z, eval );
      gaueval_spherical_angular_1_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 2:
      } else if( l == 2 ) {
        gaueval_spherical_angular_2( bf, x, y, z, eval );
      gaueval_spherical_angular_2_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 3:
      } else if( l == 3 ) {
        gaueval_spherical_angular_3( bf, x, y, z, eval );
      gaueval_spherical_angular_3_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 4:
      } else if( l == 4 ) {
        gaueval_spherical_angular_4( bf, x, y, z, eval );
      gaueval_spherical_angular_4_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 5:
      } else if( l == 5 ) {
        gaueval_spherical_angular_5( bf, x, y, z, eval );
      gaueval_spherical_angular_5_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  case 6:
      } else if( l == 6 ) {
        gaueval_spherical_angular_6( bf, x, y, z, eval );
      gaueval_spherical_angular_6_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );
//    break;

//  default: 
    } else {
      assert( false && "L < L_MAX" );
//    break;
    }

//} // switch(l)

} // gaueval_spherical_angular_deriv1



} // namespace gpgaueval
