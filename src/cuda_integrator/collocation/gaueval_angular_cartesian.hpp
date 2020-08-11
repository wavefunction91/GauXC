#pragma once
#include "gaueval_device_constants.hpp"
#include <cassert>

#ifndef GPGAUEVAL_INLINE
#  define GPGAUEVAL_INLINE __inline__
#endif

namespace gpgaueval {

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_0(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf;

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_0_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_1(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_1_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_2(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_2_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_3(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_3_deriv1(
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

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_4(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*x*x*x;
  eval[1] = bf*x*x*x*y;
  eval[2] = bf*x*x*x*z;
  eval[3] = bf*x*x*y*y;
  eval[4] = bf*x*x*y*z;
  eval[5] = bf*x*x*z*z;
  eval[6] = bf*x*y*y*y;
  eval[7] = bf*x*y*y*z;
  eval[8] = bf*x*y*z*z;
  eval[9] = bf*x*z*z*z;
  eval[10] = bf*y*y*y*y;
  eval[11] = bf*y*y*y*z;
  eval[12] = bf*y*y*z*z;
  eval[13] = bf*y*z*z*z;
  eval[14] = bf*z*z*z*z;

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_4_deriv1(
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

  eval_x[0] = x*x*x*(4*bf + bf_x*x);
  eval_x[1] = x*x*y*(3*bf + bf_x*x);
  eval_x[2] = x*x*z*(3*bf + bf_x*x);
  eval_x[3] = x*y*y*(2*bf + bf_x*x);
  eval_x[4] = x*y*z*(2*bf + bf_x*x);
  eval_x[5] = x*z*z*(2*bf + bf_x*x);
  eval_x[6] = y*y*y*(bf + bf_x*x);
  eval_x[7] = y*y*z*(bf + bf_x*x);
  eval_x[8] = y*z*z*(bf + bf_x*x);
  eval_x[9] = z*z*z*(bf + bf_x*x);
  eval_x[10] = bf_x*y*y*y*y;
  eval_x[11] = bf_x*y*y*y*z;
  eval_x[12] = bf_x*y*y*z*z;
  eval_x[13] = bf_x*y*z*z*z;
  eval_x[14] = bf_x*z*z*z*z;

  eval_y[0] = bf_y*x*x*x*x;
  eval_y[1] = x*x*x*(bf + bf_y*y);
  eval_y[2] = bf_y*x*x*x*z;
  eval_y[3] = x*x*y*(2*bf + bf_y*y);
  eval_y[4] = x*x*z*(bf + bf_y*y);
  eval_y[5] = bf_y*x*x*z*z;
  eval_y[6] = x*y*y*(3*bf + bf_y*y);
  eval_y[7] = x*y*z*(2*bf + bf_y*y);
  eval_y[8] = x*z*z*(bf + bf_y*y);
  eval_y[9] = bf_y*x*z*z*z;
  eval_y[10] = y*y*y*(4*bf + bf_y*y);
  eval_y[11] = y*y*z*(3*bf + bf_y*y);
  eval_y[12] = y*z*z*(2*bf + bf_y*y);
  eval_y[13] = z*z*z*(bf + bf_y*y);
  eval_y[14] = bf_y*z*z*z*z;

  eval_z[0] = bf_z*x*x*x*x;
  eval_z[1] = bf_z*x*x*x*y;
  eval_z[2] = x*x*x*(bf + bf_z*z);
  eval_z[3] = bf_z*x*x*y*y;
  eval_z[4] = x*x*y*(bf + bf_z*z);
  eval_z[5] = x*x*z*(2*bf + bf_z*z);
  eval_z[6] = bf_z*x*y*y*y;
  eval_z[7] = x*y*y*(bf + bf_z*z);
  eval_z[8] = x*y*z*(2*bf + bf_z*z);
  eval_z[9] = x*z*z*(3*bf + bf_z*z);
  eval_z[10] = bf_z*y*y*y*y;
  eval_z[11] = y*y*y*(bf + bf_z*z);
  eval_z[12] = y*y*z*(2*bf + bf_z*z);
  eval_z[13] = y*z*z*(3*bf + bf_z*z);
  eval_z[14] = z*z*z*(4*bf + bf_z*z);

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_5(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*x*x*x*x;
  eval[1] = bf*x*x*x*x*y;
  eval[2] = bf*x*x*x*x*z;
  eval[3] = bf*x*x*x*y*y;
  eval[4] = bf*x*x*x*y*z;
  eval[5] = bf*x*x*x*z*z;
  eval[6] = bf*x*x*y*y*y;
  eval[7] = bf*x*x*y*y*z;
  eval[8] = bf*x*x*y*z*z;
  eval[9] = bf*x*x*z*z*z;
  eval[10] = bf*x*y*y*y*y;
  eval[11] = bf*x*y*y*y*z;
  eval[12] = bf*x*y*y*z*z;
  eval[13] = bf*x*y*z*z*z;
  eval[14] = bf*x*z*z*z*z;
  eval[15] = bf*y*y*y*y*y;
  eval[16] = bf*y*y*y*y*z;
  eval[17] = bf*y*y*y*z*z;
  eval[18] = bf*y*y*z*z*z;
  eval[19] = bf*y*z*z*z*z;
  eval[20] = bf*z*z*z*z*z;

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_5_deriv1(
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

  eval_x[0] = x*x*x*x*(5*bf + bf_x*x);
  eval_x[1] = x*x*x*y*(4*bf + bf_x*x);
  eval_x[2] = x*x*x*z*(4*bf + bf_x*x);
  eval_x[3] = x*x*y*y*(3*bf + bf_x*x);
  eval_x[4] = x*x*y*z*(3*bf + bf_x*x);
  eval_x[5] = x*x*z*z*(3*bf + bf_x*x);
  eval_x[6] = x*y*y*y*(2*bf + bf_x*x);
  eval_x[7] = x*y*y*z*(2*bf + bf_x*x);
  eval_x[8] = x*y*z*z*(2*bf + bf_x*x);
  eval_x[9] = x*z*z*z*(2*bf + bf_x*x);
  eval_x[10] = y*y*y*y*(bf + bf_x*x);
  eval_x[11] = y*y*y*z*(bf + bf_x*x);
  eval_x[12] = y*y*z*z*(bf + bf_x*x);
  eval_x[13] = y*z*z*z*(bf + bf_x*x);
  eval_x[14] = z*z*z*z*(bf + bf_x*x);
  eval_x[15] = bf_x*y*y*y*y*y;
  eval_x[16] = bf_x*y*y*y*y*z;
  eval_x[17] = bf_x*y*y*y*z*z;
  eval_x[18] = bf_x*y*y*z*z*z;
  eval_x[19] = bf_x*y*z*z*z*z;
  eval_x[20] = bf_x*z*z*z*z*z;

  eval_y[0] = bf_y*x*x*x*x*x;
  eval_y[1] = x*x*x*x*(bf + bf_y*y);
  eval_y[2] = bf_y*x*x*x*x*z;
  eval_y[3] = x*x*x*y*(2*bf + bf_y*y);
  eval_y[4] = x*x*x*z*(bf + bf_y*y);
  eval_y[5] = bf_y*x*x*x*z*z;
  eval_y[6] = x*x*y*y*(3*bf + bf_y*y);
  eval_y[7] = x*x*y*z*(2*bf + bf_y*y);
  eval_y[8] = x*x*z*z*(bf + bf_y*y);
  eval_y[9] = bf_y*x*x*z*z*z;
  eval_y[10] = x*y*y*y*(4*bf + bf_y*y);
  eval_y[11] = x*y*y*z*(3*bf + bf_y*y);
  eval_y[12] = x*y*z*z*(2*bf + bf_y*y);
  eval_y[13] = x*z*z*z*(bf + bf_y*y);
  eval_y[14] = bf_y*x*z*z*z*z;
  eval_y[15] = y*y*y*y*(5*bf + bf_y*y);
  eval_y[16] = y*y*y*z*(4*bf + bf_y*y);
  eval_y[17] = y*y*z*z*(3*bf + bf_y*y);
  eval_y[18] = y*z*z*z*(2*bf + bf_y*y);
  eval_y[19] = z*z*z*z*(bf + bf_y*y);
  eval_y[20] = bf_y*z*z*z*z*z;

  eval_z[0] = bf_z*x*x*x*x*x;
  eval_z[1] = bf_z*x*x*x*x*y;
  eval_z[2] = x*x*x*x*(bf + bf_z*z);
  eval_z[3] = bf_z*x*x*x*y*y;
  eval_z[4] = x*x*x*y*(bf + bf_z*z);
  eval_z[5] = x*x*x*z*(2*bf + bf_z*z);
  eval_z[6] = bf_z*x*x*y*y*y;
  eval_z[7] = x*x*y*y*(bf + bf_z*z);
  eval_z[8] = x*x*y*z*(2*bf + bf_z*z);
  eval_z[9] = x*x*z*z*(3*bf + bf_z*z);
  eval_z[10] = bf_z*x*y*y*y*y;
  eval_z[11] = x*y*y*y*(bf + bf_z*z);
  eval_z[12] = x*y*y*z*(2*bf + bf_z*z);
  eval_z[13] = x*y*z*z*(3*bf + bf_z*z);
  eval_z[14] = x*z*z*z*(4*bf + bf_z*z);
  eval_z[15] = bf_z*y*y*y*y*y;
  eval_z[16] = y*y*y*y*(bf + bf_z*z);
  eval_z[17] = y*y*y*z*(2*bf + bf_z*z);
  eval_z[18] = y*y*z*z*(3*bf + bf_z*z);
  eval_z[19] = y*z*z*z*(4*bf + bf_z*z);
  eval_z[20] = z*z*z*z*(5*bf + bf_z*z);

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_6(
  const double bf,
  const double x,
  const double y,
  const double z,
  double*      eval
) {

  eval[0] = bf*x*x*x*x*x*x;
  eval[1] = bf*x*x*x*x*x*y;
  eval[2] = bf*x*x*x*x*x*z;
  eval[3] = bf*x*x*x*x*y*y;
  eval[4] = bf*x*x*x*x*y*z;
  eval[5] = bf*x*x*x*x*z*z;
  eval[6] = bf*x*x*x*y*y*y;
  eval[7] = bf*x*x*x*y*y*z;
  eval[8] = bf*x*x*x*y*z*z;
  eval[9] = bf*x*x*x*z*z*z;
  eval[10] = bf*x*x*y*y*y*y;
  eval[11] = bf*x*x*y*y*y*z;
  eval[12] = bf*x*x*y*y*z*z;
  eval[13] = bf*x*x*y*z*z*z;
  eval[14] = bf*x*x*z*z*z*z;
  eval[15] = bf*x*y*y*y*y*y;
  eval[16] = bf*x*y*y*y*y*z;
  eval[17] = bf*x*y*y*y*z*z;
  eval[18] = bf*x*y*y*z*z*z;
  eval[19] = bf*x*y*z*z*z*z;
  eval[20] = bf*x*z*z*z*z*z;
  eval[21] = bf*y*y*y*y*y*y;
  eval[22] = bf*y*y*y*y*y*z;
  eval[23] = bf*y*y*y*y*z*z;
  eval[24] = bf*y*y*y*z*z*z;
  eval[25] = bf*y*y*z*z*z*z;
  eval[26] = bf*y*z*z*z*z*z;
  eval[27] = bf*z*z*z*z*z*z;

}

GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_6_deriv1(
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

  eval_x[0] = x*x*x*x*x*(6*bf + bf_x*x);
  eval_x[1] = x*x*x*x*y*(5*bf + bf_x*x);
  eval_x[2] = x*x*x*x*z*(5*bf + bf_x*x);
  eval_x[3] = x*x*x*y*y*(4*bf + bf_x*x);
  eval_x[4] = x*x*x*y*z*(4*bf + bf_x*x);
  eval_x[5] = x*x*x*z*z*(4*bf + bf_x*x);
  eval_x[6] = x*x*y*y*y*(3*bf + bf_x*x);
  eval_x[7] = x*x*y*y*z*(3*bf + bf_x*x);
  eval_x[8] = x*x*y*z*z*(3*bf + bf_x*x);
  eval_x[9] = x*x*z*z*z*(3*bf + bf_x*x);
  eval_x[10] = x*y*y*y*y*(2*bf + bf_x*x);
  eval_x[11] = x*y*y*y*z*(2*bf + bf_x*x);
  eval_x[12] = x*y*y*z*z*(2*bf + bf_x*x);
  eval_x[13] = x*y*z*z*z*(2*bf + bf_x*x);
  eval_x[14] = x*z*z*z*z*(2*bf + bf_x*x);
  eval_x[15] = y*y*y*y*y*(bf + bf_x*x);
  eval_x[16] = y*y*y*y*z*(bf + bf_x*x);
  eval_x[17] = y*y*y*z*z*(bf + bf_x*x);
  eval_x[18] = y*y*z*z*z*(bf + bf_x*x);
  eval_x[19] = y*z*z*z*z*(bf + bf_x*x);
  eval_x[20] = z*z*z*z*z*(bf + bf_x*x);
  eval_x[21] = bf_x*y*y*y*y*y*y;
  eval_x[22] = bf_x*y*y*y*y*y*z;
  eval_x[23] = bf_x*y*y*y*y*z*z;
  eval_x[24] = bf_x*y*y*y*z*z*z;
  eval_x[25] = bf_x*y*y*z*z*z*z;
  eval_x[26] = bf_x*y*z*z*z*z*z;
  eval_x[27] = bf_x*z*z*z*z*z*z;

  eval_y[0] = bf_y*x*x*x*x*x*x;
  eval_y[1] = x*x*x*x*x*(bf + bf_y*y);
  eval_y[2] = bf_y*x*x*x*x*x*z;
  eval_y[3] = x*x*x*x*y*(2*bf + bf_y*y);
  eval_y[4] = x*x*x*x*z*(bf + bf_y*y);
  eval_y[5] = bf_y*x*x*x*x*z*z;
  eval_y[6] = x*x*x*y*y*(3*bf + bf_y*y);
  eval_y[7] = x*x*x*y*z*(2*bf + bf_y*y);
  eval_y[8] = x*x*x*z*z*(bf + bf_y*y);
  eval_y[9] = bf_y*x*x*x*z*z*z;
  eval_y[10] = x*x*y*y*y*(4*bf + bf_y*y);
  eval_y[11] = x*x*y*y*z*(3*bf + bf_y*y);
  eval_y[12] = x*x*y*z*z*(2*bf + bf_y*y);
  eval_y[13] = x*x*z*z*z*(bf + bf_y*y);
  eval_y[14] = bf_y*x*x*z*z*z*z;
  eval_y[15] = x*y*y*y*y*(5*bf + bf_y*y);
  eval_y[16] = x*y*y*y*z*(4*bf + bf_y*y);
  eval_y[17] = x*y*y*z*z*(3*bf + bf_y*y);
  eval_y[18] = x*y*z*z*z*(2*bf + bf_y*y);
  eval_y[19] = x*z*z*z*z*(bf + bf_y*y);
  eval_y[20] = bf_y*x*z*z*z*z*z;
  eval_y[21] = y*y*y*y*y*(6*bf + bf_y*y);
  eval_y[22] = y*y*y*y*z*(5*bf + bf_y*y);
  eval_y[23] = y*y*y*z*z*(4*bf + bf_y*y);
  eval_y[24] = y*y*z*z*z*(3*bf + bf_y*y);
  eval_y[25] = y*z*z*z*z*(2*bf + bf_y*y);
  eval_y[26] = z*z*z*z*z*(bf + bf_y*y);
  eval_y[27] = bf_y*z*z*z*z*z*z;

  eval_z[0] = bf_z*x*x*x*x*x*x;
  eval_z[1] = bf_z*x*x*x*x*x*y;
  eval_z[2] = x*x*x*x*x*(bf + bf_z*z);
  eval_z[3] = bf_z*x*x*x*x*y*y;
  eval_z[4] = x*x*x*x*y*(bf + bf_z*z);
  eval_z[5] = x*x*x*x*z*(2*bf + bf_z*z);
  eval_z[6] = bf_z*x*x*x*y*y*y;
  eval_z[7] = x*x*x*y*y*(bf + bf_z*z);
  eval_z[8] = x*x*x*y*z*(2*bf + bf_z*z);
  eval_z[9] = x*x*x*z*z*(3*bf + bf_z*z);
  eval_z[10] = bf_z*x*x*y*y*y*y;
  eval_z[11] = x*x*y*y*y*(bf + bf_z*z);
  eval_z[12] = x*x*y*y*z*(2*bf + bf_z*z);
  eval_z[13] = x*x*y*z*z*(3*bf + bf_z*z);
  eval_z[14] = x*x*z*z*z*(4*bf + bf_z*z);
  eval_z[15] = bf_z*x*y*y*y*y*y;
  eval_z[16] = x*y*y*y*y*(bf + bf_z*z);
  eval_z[17] = x*y*y*y*z*(2*bf + bf_z*z);
  eval_z[18] = x*y*y*z*z*(3*bf + bf_z*z);
  eval_z[19] = x*y*z*z*z*(4*bf + bf_z*z);
  eval_z[20] = x*z*z*z*z*(5*bf + bf_z*z);
  eval_z[21] = bf_z*y*y*y*y*y*y;
  eval_z[22] = y*y*y*y*y*(bf + bf_z*z);
  eval_z[23] = y*y*y*y*z*(2*bf + bf_z*z);
  eval_z[24] = y*y*y*z*z*(3*bf + bf_z*z);
  eval_z[25] = y*y*z*z*z*(4*bf + bf_z*z);
  eval_z[26] = y*z*z*z*z*(5*bf + bf_z*z);
  eval_z[27] = z*z*z*z*z*(6*bf + bf_z*z);

}


GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular(
  const int64_t l,
  const double  bf,
  const double  x,
  const double  y,
  const double  z,
  double*       eval
) {

      if( l == 0 ) {
        gaueval_cartesian_angular_0( bf, x, y, z, eval );

      } else if( l == 1 ) {
        gaueval_cartesian_angular_1( bf, x, y, z, eval );

      } else if( l == 2 ) {
        gaueval_cartesian_angular_2( bf, x, y, z, eval );

      } else if( l == 3 ) {
        gaueval_cartesian_angular_3( bf, x, y, z, eval );

      } else if( l == 4 ) {
        gaueval_cartesian_angular_4( bf, x, y, z, eval );

      } else if( l == 5 ) {
        gaueval_cartesian_angular_5( bf, x, y, z, eval );

      } else if( l == 6 ) {
        gaueval_cartesian_angular_6( bf, x, y, z, eval );

    } else {
      assert( false && "L < L_MAX" );
    }

} // gaueval_cartesian_angular


GPGAUEVAL_INLINE __device__ void gaueval_cartesian_angular_deriv1(
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
        gaueval_cartesian_angular_0( bf, x, y, z, eval );
      gaueval_cartesian_angular_0_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 1 ) {
        gaueval_cartesian_angular_1( bf, x, y, z, eval );
      gaueval_cartesian_angular_1_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 2 ) {
        gaueval_cartesian_angular_2( bf, x, y, z, eval );
      gaueval_cartesian_angular_2_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 3 ) {
        gaueval_cartesian_angular_3( bf, x, y, z, eval );
      gaueval_cartesian_angular_3_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 4 ) {
        gaueval_cartesian_angular_4( bf, x, y, z, eval );
      gaueval_cartesian_angular_4_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 5 ) {
        gaueval_cartesian_angular_5( bf, x, y, z, eval );
      gaueval_cartesian_angular_5_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

      } else if( l == 6 ) {
        gaueval_cartesian_angular_6( bf, x, y, z, eval );
      gaueval_cartesian_angular_6_deriv1( bf, bf_x, bf_y, bf_z, x, y, z, eval_x, eval_y, eval_z );

    } else {
      assert( false && "L < L_MAX" );
    }

} // gaueval_cartesian_angular_deriv1



} // namespace gpgaueval
