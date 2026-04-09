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
#include "device/common/onedft_exc_grad.hpp"
#include <gauxc/util/div_ceil.hpp>
#include "device_specific/cuda_util.hpp"
#include "device_specific/cuda_device_constants.hpp"

namespace GauXC {

// Transform OneDFT per-component Vxc to the standard format expected by
// inc_exc_grad_gga/mgga kernels. In OneDFT:
//   gamma_pp  = dden_x_grad_a  (per-direction derivative, alpha)
//   vgamma_pp = dden_x_grad_b  (per-direction derivative, beta)
//   gamma_pm  = dden_y_grad_a
//   vgamma_pm = dden_y_grad_b
//   gamma_mm  = dden_z_grad_a
//   vgamma_mm = dden_z_grad_b
//
// The standard kernel reads dden_sx = (a+b)/2 (s-spin), dden_zx = (a-b)/2 (z-spin)
// and vgamma_pp/pm/mm as scalar coupling coefficients.
// Setting vgamma_pp=1, vgamma_pm=0, vgamma_mm=1 makes the standard kernel
// reproduce the OneDFT gradient formula exactly.
__global__ void transform_onedft_vxc_for_grad_kernel(
    uint32_t ntasks,
    XCDeviceTask* __restrict__ tasks_device ) {

  const int batch_idx = blockIdx.z;
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[batch_idx];
  const auto npts = task.npts;

  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if( tid >= npts ) return;

  // Read per-direction OneDFT derivatives (alpha/beta)
  const double dx_a = task.gamma_pp[tid];
  const double dx_b = task.vgamma_pp[tid];
  const double dy_a = task.gamma_pm[tid];
  const double dy_b = task.vgamma_pm[tid];
  const double dz_a = task.gamma_mm[tid];
  const double dz_b = task.vgamma_mm[tid];

  // Convert to total (s) and magnetization (z) form
  task.dden_sx[tid] = 0.5 * (dx_a + dx_b);
  task.dden_sy[tid] = 0.5 * (dy_a + dy_b);
  task.dden_sz[tid] = 0.5 * (dz_a + dz_b);
  task.dden_zx[tid] = 0.5 * (dx_a - dx_b);
  task.dden_zy[tid] = 0.5 * (dy_a - dy_b);
  task.dden_zz[tid] = 0.5 * (dz_a - dz_b);

  // Set vgamma coefficients so the standard kernel reproduces the OneDFT formula
  task.vgamma_pp[tid] = 1.0;
  task.vgamma_pm[tid] = 0.0;
  task.vgamma_mm[tid] = 1.0;
}

void transform_onedft_vxc_for_grad(
    size_t ntasks,
    int32_t max_npts,
    XCDeviceTask* tasks_device,
    device_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();

  dim3 threads(256);
  dim3 blocks( util::div_ceil((uint32_t)max_npts, threads.x), 1, ntasks );

  transform_onedft_vxc_for_grad_kernel<<<blocks, threads, 0, stream>>>(
    ntasks, tasks_device );
}

} // namespace GauXC
