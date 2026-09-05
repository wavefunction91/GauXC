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
#include "sycl_aos_scheme1.hpp"
#include "device/sycl/sycl_backend.hpp"
#include "kernels/grid_to_center.hpp"
#include "kernels/sycl_ssf_2d.hpp"
#include "kernels/sycl_ssf_1d.hpp"
#include "kernels/sycl_launch.hpp"
#include "sycl_aos_scheme1_weights.hpp"
 
namespace GauXC {

void sycl_aos_scheme1_weights_wrapper( int32_t npts, int32_t natoms,
  const double* points_x, const double* points_y, const double* points_z, 
  const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, double* weights, ::sycl::queue& stream ) {

  // Compute distances from grid to atomic centers
  compute_grid_to_center_dist( npts, natoms, coords, points_x, points_y, points_z, 
   dist, lddist, stream );

  // Modify weights. The 2D kernel is retained (and exercised below) because,
  // unlike the CUDA backend, PVC has no XeCore-count-driven persistent-block
  // launch to tune against; the 1D kernel is the portable default
  partition_weights_ssf_1d( npts, natoms, RAB, ldRAB, coords, dist, lddist,
    iparent, dist_nearest, weights, stream );

}


void sycl_aos_scheme1_weight_1st_deriv_wrapper(
  int32_t npts, int32_t natoms,
  const double* points_x, const double* points_y, const double* points_z,
  const double* RAB, int32_t ldRAB, const double* coords, 
  double* dist, int32_t lddist, const int32_t* iparent,
  const double* dist_nearest, const double* w_times_f,
  double* exc_grad_w, ::sycl::queue& stream ){

  // Compute distances from grid to atomic centers
  compute_grid_to_center_dist( npts, natoms, coords, points_x, points_y, points_z, 
   dist, lddist, stream );

  eval_weight_1st_deriv_contracted_ssf_1d( npts, natoms, RAB, ldRAB, coords,
    points_x, points_y, points_z, dist, lddist, iparent, dist_nearest,
    w_times_f, exc_grad_w, stream );

}





}
