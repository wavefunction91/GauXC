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

#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_HAS_CUDA

#include "../host/vv10_nlc.hpp"
#include "device/cuda/kernels/vv10.hpp"
#include "device_specific/cuda_util.hpp"

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace GauXC::detail::vv10 {

inline Parameters make_cuda_parameters( const IntegratorSettingsNLC& settings ) {
  return Parameters{ settings.vv10_b, settings.vv10_c, settings.vv10_tol,
    settings.math_mode == NLCMathMode::FloatPair ? 1 : 0 };
}

inline void eval_exc_vxc_cuda( cudaStream_t stream, const IntegratorSettingsNLC& settings,
                               GridView grid, CorrectionsView corr ) {
  if( grid.npts == 0 ) return;

  std::vector<std::size_t> active_indices;
  active_indices.reserve(grid.npts);
  for( std::size_t i = 0; i < grid.npts; ++i ) {
    if( grid.rho[i] >= settings.vv10_tol ) active_indices.emplace_back(i);
  }
  if( active_indices.empty() ) return;

  const auto active_npts = active_indices.size();
  const auto npts = static_cast<int32_t>( active_npts );
  const bool all_active = active_npts == grid.npts;
  std::vector<double> coords, rho, gamma, rho_weight(active_npts);
  const double* coords_host = grid.coords;
  const double* rho_host = grid.rho;
  const double* gamma_host = grid.gamma;
  if( all_active ) {
    for( std::size_t i = 0; i < active_npts; ++i ) {
      rho_weight[i] = grid.rho[i] * grid.weights[i];
    }
  } else {
    coords.resize(3 * active_npts);
    rho.resize(active_npts);
    gamma.resize(active_npts);
    for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
      const auto i = active_indices[active_i];
      coords[3*active_i + 0] = grid.coords[3*i + 0];
      coords[3*active_i + 1] = grid.coords[3*i + 1];
      coords[3*active_i + 2] = grid.coords[3*i + 2];
      rho[active_i] = grid.rho[i];
      gamma[active_i] = grid.gamma[i];
      rho_weight[active_i] = grid.rho[i] * grid.weights[i];
    }
    coords_host = coords.data();
    rho_host = rho.data();
    gamma_host = gamma.data();
  }

  auto* coords_device     = util::cuda_malloc<double>( 3 * active_npts );
  auto* rho_device        = util::cuda_malloc<double>( active_npts );
  auto* gamma_device      = util::cuda_malloc<double>( active_npts );
  auto* rho_weight_device = util::cuda_malloc<double>( active_npts );
  auto* omega_device      = util::cuda_malloc<double>( active_npts );
  auto* kappa_device      = util::cuda_malloc<double>( active_npts );
  auto* eps_device        = util::cuda_malloc<double>( active_npts );
  auto* vrho_device       = util::cuda_malloc<double>( active_npts );
  auto* vgamma_device     = util::cuda_malloc<double>( active_npts );

  util::cuda_copy( 3 * active_npts, coords_device, coords_host );
  util::cuda_copy( active_npts, rho_device, rho_host );
  util::cuda_copy( active_npts, gamma_device, gamma_host );
  util::cuda_copy( active_npts, rho_weight_device, rho_weight.data() );

  const auto params = make_cuda_parameters(settings);
  eval_omega_kappa_cuda( stream, npts, rho_device, gamma_device, params,
    omega_device, kappa_device, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr );
  eval_exc_vxc_fock_cuda( stream, npts, coords_device, rho_device, gamma_device,
    omega_device, kappa_device, rho_weight_device, params, eps_device,
    vrho_device, vgamma_device );
  util::cuda_device_sync();

  std::vector<double> eps(active_npts), vrho(active_npts), vgamma(active_npts);
  util::cuda_copy( active_npts, eps.data(), eps_device );
  util::cuda_copy( active_npts, vrho.data(), vrho_device );
  util::cuda_copy( active_npts, vgamma.data(), vgamma_device );

  if( all_active ) {
    for( std::size_t i = 0; i < active_npts; ++i ) {
      corr.eps[i]    += eps[i];
      corr.vrho[i]   += vrho[i];
      corr.vgamma[i] += vgamma[i];
    }
  } else {
    for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
      const auto i = active_indices[active_i];
      corr.eps[i]    += eps[active_i];
      corr.vrho[i]   += vrho[active_i];
      corr.vgamma[i] += vgamma[active_i];
    }
  }

  util::cuda_free( coords_device, rho_device, gamma_device, rho_weight_device,
    omega_device, kappa_device, eps_device, vrho_device, vgamma_device );
}

inline void eval_fxc_rks_cuda( cudaStream_t stream, const IntegratorSettingsNLC& settings,
                               ResponseGridView grid, TrialView trial,
                               ResponseCorrectionsView corr ) {
  if( grid.npts == 0 ) return;

  std::vector<std::size_t> active_indices;
  active_indices.reserve(grid.npts);
  for( std::size_t i = 0; i < grid.npts; ++i ) {
    if( grid.rho[i] >= settings.vv10_tol ) active_indices.emplace_back(i);
  }
  if( active_indices.empty() ) return;

  const auto active_npts = active_indices.size();
  const auto active_npts_i32 = static_cast<int32_t>(active_npts);
  std::vector<double> coords(3 * active_npts), weights(active_npts), rho(active_npts), gamma(active_npts);
  std::vector<double> rho_t(active_npts), gamma_t(active_npts);

  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    coords[3*active_i + 0] = grid.coords[3*i + 0];
    coords[3*active_i + 1] = grid.coords[3*i + 1];
    coords[3*active_i + 2] = grid.coords[3*i + 2];
    weights[active_i] = grid.weights[i];
    rho[active_i] = grid.rho[i];
    gamma[active_i] = grid.gamma[i];
    rho_t[active_i] = trial.rho[i];
    gamma_t[active_i] = 2.0 * (grid.grad_x[i] * trial.grad_x[i]
      + grid.grad_y[i] * trial.grad_y[i] + grid.grad_z[i] * trial.grad_z[i]);
  }

  auto* coords_device = util::cuda_malloc<double>(3 * active_npts);
  auto* weights_device = util::cuda_malloc<double>(active_npts);
  auto* rho_device = util::cuda_malloc<double>(active_npts);
  auto* gamma_device = util::cuda_malloc<double>(active_npts);
  auto* omega_device = util::cuda_malloc<double>(active_npts);
  auto* kappa_device = util::cuda_malloc<double>(active_npts);
  auto* domega_drho_device = util::cuda_malloc<double>(active_npts);
  auto* domega_dgamma_device = util::cuda_malloc<double>(active_npts);
  auto* d2omega_drho2_device = util::cuda_malloc<double>(active_npts);
  auto* d2omega_dgamma2_device = util::cuda_malloc<double>(active_npts);
  auto* d2omega_drho_dgamma_device = util::cuda_malloc<double>(active_npts);
  auto* dkappa_drho_device = util::cuda_malloc<double>(active_npts);
  auto* d2kappa_drho2_device = util::cuda_malloc<double>(active_npts);
  auto* u_device = util::cuda_malloc<double>(active_npts);
  auto* w_device = util::cuda_malloc<double>(active_npts);
  auto* a_device = util::cuda_malloc<double>(active_npts);
  auto* b_device = util::cuda_malloc<double>(active_npts);
  auto* c_device = util::cuda_malloc<double>(active_npts);
  auto* e_device = util::cuda_malloc<double>(active_npts);
  auto* rho_t_device = util::cuda_malloc<double>(active_npts);
  auto* gamma_t_device = util::cuda_malloc<double>(active_npts);
  auto* f_rho_t_device = util::cuda_malloc<double>(active_npts);
  auto* f_gamma_t_device = util::cuda_malloc<double>(active_npts);

  util::cuda_copy(3 * active_npts, coords_device, coords.data());
  util::cuda_copy(active_npts, weights_device, weights.data());
  util::cuda_copy(active_npts, rho_device, rho.data());
  util::cuda_copy(active_npts, gamma_device, gamma.data());
  util::cuda_copy(active_npts, rho_t_device, rho_t.data());
  util::cuda_copy(active_npts, gamma_t_device, gamma_t.data());

  const auto params = make_cuda_parameters(settings);
  eval_omega_kappa_cuda(stream, active_npts_i32, rho_device, gamma_device, params,
    omega_device, kappa_device, domega_drho_device, domega_dgamma_device,
    d2omega_drho2_device, d2omega_dgamma2_device, d2omega_drho_dgamma_device,
    dkappa_drho_device, d2kappa_drho2_device);
  eval_hessian_intermediates_cuda(stream, active_npts_i32, coords_device,
    weights_device, rho_device, omega_device, kappa_device, u_device, w_device,
    a_device, b_device, c_device, e_device);
  eval_hessian_contraction_cuda(stream, active_npts_i32, 1, coords_device,
    weights_device, rho_device, omega_device, kappa_device, u_device, w_device,
    a_device, b_device, c_device, domega_drho_device, domega_dgamma_device,
    dkappa_drho_device, d2omega_drho2_device, d2omega_dgamma2_device,
    d2omega_drho_dgamma_device, d2kappa_drho2_device, rho_t_device,
    gamma_t_device, f_rho_t_device, f_gamma_t_device);
  util::cuda_device_sync();

  std::vector<double> w(active_npts), domega_dgamma(active_npts);
  std::vector<double> f_rho_t(active_npts), f_gamma_t(active_npts);
  util::cuda_copy(active_npts, w.data(), w_device);
  util::cuda_copy(active_npts, domega_dgamma.data(), domega_dgamma_device);
  util::cuda_copy(active_npts, f_rho_t.data(), f_rho_t_device);
  util::cuda_copy(active_npts, f_gamma_t.data(), f_gamma_t_device);

  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    const double f_gamma = grid.rho[i] * domega_dgamma[active_i] * w[active_i];
    corr.A[i] += f_rho_t[active_i];
    corr.B[3*i]   += 2.0 * (grid.grad_x[i] * f_gamma_t[active_i] + trial.grad_x[i] * f_gamma);
    corr.B[3*i+1] += 2.0 * (grid.grad_y[i] * f_gamma_t[active_i] + trial.grad_y[i] * f_gamma);
    corr.B[3*i+2] += 2.0 * (grid.grad_z[i] * f_gamma_t[active_i] + trial.grad_z[i] * f_gamma);
  }

  util::cuda_free(coords_device, weights_device, rho_device, gamma_device,
    omega_device, kappa_device, domega_drho_device, domega_dgamma_device,
    d2omega_drho2_device, d2omega_dgamma2_device, d2omega_drho_dgamma_device,
    dkappa_drho_device, d2kappa_drho2_device, u_device, w_device, a_device,
    b_device, c_device, e_device, rho_t_device, gamma_t_device,
    f_rho_t_device, f_gamma_t_device);
}

inline void eval_grid_gradient_cuda( cudaStream_t stream, const IntegratorSettingsNLC& settings,
                                     GridView grid, GridGradientView grad ) {
  if( grid.npts == 0 ) return;

  std::vector<std::size_t> active_indices;
  active_indices.reserve(grid.npts);
  for( std::size_t i = 0; i < grid.npts; ++i ) {
    if( grid.rho[i] >= settings.vv10_tol ) active_indices.emplace_back(i);
  }
  if( active_indices.empty() ) return;

  const auto active_npts = active_indices.size();
  const auto active_npts_i32 = static_cast<int32_t>(active_npts);
  std::vector<double> coords(3 * active_npts), rho(active_npts), gamma(active_npts), rho_weight(active_npts);
  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    coords[3*active_i + 0] = grid.coords[3*i + 0];
    coords[3*active_i + 1] = grid.coords[3*i + 1];
    coords[3*active_i + 2] = grid.coords[3*i + 2];
    rho[active_i] = grid.rho[i];
    gamma[active_i] = grid.gamma[i];
    rho_weight[active_i] = grid.rho[i] * grid.weights[i];
  }

  auto* coords_device = util::cuda_malloc<double>(3 * active_npts);
  auto* rho_device = util::cuda_malloc<double>(active_npts);
  auto* gamma_device = util::cuda_malloc<double>(active_npts);
  auto* rho_weight_device = util::cuda_malloc<double>(active_npts);
  auto* omega_device = util::cuda_malloc<double>(active_npts);
  auto* kappa_device = util::cuda_malloc<double>(active_npts);
  auto* gradient_device = util::cuda_malloc<double>(3 * active_npts);

  util::cuda_copy(3 * active_npts, coords_device, coords.data());
  util::cuda_copy(active_npts, rho_device, rho.data());
  util::cuda_copy(active_npts, gamma_device, gamma.data());
  util::cuda_copy(active_npts, rho_weight_device, rho_weight.data());

  const auto params = make_cuda_parameters(settings);
  eval_omega_kappa_cuda(stream, active_npts_i32, rho_device, gamma_device, params,
    omega_device, kappa_device, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  eval_grid_gradient_cuda(stream, active_npts_i32, coords_device, omega_device,
    kappa_device, active_npts_i32, coords_device, rho_weight_device,
    omega_device, kappa_device, gradient_device);
  util::cuda_device_sync();

  std::vector<double> gradient(3 * active_npts);
  util::cuda_copy(3 * active_npts, gradient.data(), gradient_device);
  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    grad.x[i] += gradient[3*active_i + 0];
    grad.y[i] += gradient[3*active_i + 1];
    grad.z[i] += gradient[3*active_i + 2];
  }

  util::cuda_free(coords_device, rho_device, gamma_device, rho_weight_device,
    omega_device, kappa_device, gradient_device);
}

inline void eval_grid_gradient_excluding_same_parent_cuda( cudaStream_t stream,
                                                           const IntegratorSettingsNLC& settings,
                                                           GridView grid,
                                                           const std::vector<unsigned long long>& parent,
                                                           GridGradientView grad ) {
  if( grid.npts == 0 ) return;
  if( parent.size() != grid.npts ) {
    GAUXC_GENERIC_EXCEPTION("Invalid VV10 device parent grid size");
  }

  std::vector<std::size_t> active_indices;
  active_indices.reserve(grid.npts);
  for( std::size_t i = 0; i < grid.npts; ++i ) {
    if( grid.rho[i] >= settings.vv10_tol ) active_indices.emplace_back(i);
  }
  if( active_indices.empty() ) return;

  const auto active_npts = active_indices.size();
  const auto active_npts_i32 = static_cast<int32_t>(active_npts);
  std::vector<double> coords(3 * active_npts), rho(active_npts), gamma(active_npts), rho_weight(active_npts);
  std::vector<std::int32_t> active_parent(active_npts);
  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    coords[3*active_i + 0] = grid.coords[3*i + 0];
    coords[3*active_i + 1] = grid.coords[3*i + 1];
    coords[3*active_i + 2] = grid.coords[3*i + 2];
    rho[active_i] = grid.rho[i];
    gamma[active_i] = grid.gamma[i];
    rho_weight[active_i] = grid.rho[i] * grid.weights[i];
    active_parent[active_i] = static_cast<std::int32_t>(parent[i]);
  }

  auto* coords_device = util::cuda_malloc<double>(3 * active_npts);
  auto* rho_device = util::cuda_malloc<double>(active_npts);
  auto* gamma_device = util::cuda_malloc<double>(active_npts);
  auto* rho_weight_device = util::cuda_malloc<double>(active_npts);
  auto* omega_device = util::cuda_malloc<double>(active_npts);
  auto* kappa_device = util::cuda_malloc<double>(active_npts);
  auto* parent_device = util::cuda_malloc<std::int32_t>(active_npts);
  auto* gradient_device = util::cuda_malloc<double>(3 * active_npts);

  util::cuda_copy(3 * active_npts, coords_device, coords.data());
  util::cuda_copy(active_npts, rho_device, rho.data());
  util::cuda_copy(active_npts, gamma_device, gamma.data());
  util::cuda_copy(active_npts, rho_weight_device, rho_weight.data());
  util::cuda_copy(active_npts, parent_device, active_parent.data());

  const auto params = make_cuda_parameters(settings);
  eval_omega_kappa_cuda(stream, active_npts_i32, rho_device, gamma_device, params,
    omega_device, kappa_device, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
  eval_grid_gradient_excluding_same_parent_cuda(stream, active_npts_i32,
    coords_device, omega_device, kappa_device, rho_weight_device,
    parent_device, gradient_device);
  util::cuda_device_sync();

  std::vector<double> gradient(3 * active_npts);
  util::cuda_copy(3 * active_npts, gradient.data(), gradient_device);
  for( std::size_t active_i = 0; active_i < active_npts; ++active_i ) {
    const auto i = active_indices[active_i];
    grad.x[i] += gradient[3*active_i + 0];
    grad.y[i] += gradient[3*active_i + 1];
    grad.z[i] += gradient[3*active_i + 2];
  }

  util::cuda_free(coords_device, rho_device, gamma_device, rho_weight_device,
    omega_device, kappa_device, parent_device, gradient_device);
}

} // namespace GauXC::detail::vv10

#endif