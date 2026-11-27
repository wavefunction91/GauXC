/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <sycl/sycl.hpp>

namespace GauXC {

void sycl_aos_scheme1_weights_wrapper( int32_t npts, int32_t natoms,
  const double* points_x, const double* points_y, const double* points_z,
  const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, double* weights, ::sycl::queue& stream );

void sycl_aos_scheme1_weight_1st_deriv_wrapper(
  int32_t npts, int32_t natoms,
  const double* points_x, const double* points_y, const double* points_z,
  const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, const double* w_times_f,
  double* exc_grad_w, ::sycl::queue& stream );

}
