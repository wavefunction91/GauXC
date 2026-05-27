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

#include <gauxc/exceptions.hpp>
#include <gauxc/reduction_driver.hpp>
#include <gauxc/xc_integrator_settings.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace GauXC::detail::vv10 {

struct GridView {
  std::size_t npts = 0;
  const double* coords  = nullptr;
  const double* weights = nullptr;
  const double* rho     = nullptr;
  const double* gamma   = nullptr;
};

struct CorrectionsView {
  double* eps    = nullptr;
  double* vrho   = nullptr;
  double* vgamma = nullptr;
};

struct ResponseGridView {
  std::size_t npts = 0;
  const double* coords = nullptr;
  const double* weights = nullptr;
  const double* rho = nullptr;
  const double* gamma = nullptr;
  const double* grad_x = nullptr;
  const double* grad_y = nullptr;
  const double* grad_z = nullptr;
};

struct TrialView {
  const double* rho = nullptr;
  const double* grad_x = nullptr;
  const double* grad_y = nullptr;
  const double* grad_z = nullptr;
};

struct ResponseCorrectionsView {
  double* A = nullptr;
  double* B = nullptr;
};

struct GridGradientView {
  double* x = nullptr;
  double* y = nullptr;
  double* z = nullptr;
};

inline void add_rks_correction( std::size_t i, double eps, double vrho, double vgamma,
                                double* eps_out, double* vrho_out, double* vgamma_out ) {
  eps_out[i] += eps;
  vrho_out[i] += vrho;
  vgamma_out[i] += vgamma;
}

inline void add_uks_correction( std::size_t i, double eps, double vrho, double vgamma,
                                double* eps_out, double* vrho_out, double* vgamma_out ) {
  eps_out[i] += eps;
  vrho_out[2*i]   += vrho;
  vrho_out[2*i+1] += vrho;
  vgamma_out[3*i]   += vgamma;
  vgamma_out[3*i+1] += 2.0 * vgamma;
  vgamma_out[3*i+2] += vgamma;
}

inline IntegratorSettingsNLC get_settings( const IntegratorSettingsXC& settings ) {
  IntegratorSettingsNLC vv10_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsNLC*>(&settings) ) {
    vv10_settings = *tmp;
  }
  return vv10_settings;
}

inline bool enabled( const IntegratorSettingsXC& settings ) {
  return dynamic_cast<const IntegratorSettingsNLCInternal*>(&settings) != nullptr;
}

inline double beta( const IntegratorSettingsNLC& settings ) {
  return std::pow( 3.0 / (settings.vv10_b * settings.vv10_b), 0.75 ) / 32.0;
}

inline double kappa_prefactor( const IntegratorSettingsNLC& settings ) {
  const double pi = 3.141592653589793238462643383279502884;
  return settings.vv10_b * 1.5 * pi * std::pow( 9.0 * pi, -1.0 / 6.0 );
}

inline std::vector<double> allgather_packed_grid( const ReductionDriver& reduction_driver,
                                                  const std::vector<double>& local_packed,
                                                  std::size_t stride,
                                                  std::size_t& local_point_offset ) {
  const auto local_npts = static_cast<unsigned long long>( local_packed.size() / stride );
  const auto all_npts = reduction_driver.allgather_v( &local_npts, 1 );

  local_point_offset = 0;
  const int rank = reduction_driver.comm_rank();
  for( int i = 0; i < rank; ++i ) {
    local_point_offset += static_cast<std::size_t>( all_npts[i] );
  }

  return reduction_driver.allgather_v( local_packed.data(), local_packed.size() );
}

inline void eval_exc_vxc( const IntegratorSettingsNLC& settings,
                          GridView grid, CorrectionsView corr ) {
  if( grid.npts == 0 ) return;

  const auto npts = grid.npts;
  const double pi = 3.141592653589793238462643383279502884;
  const double pi43 = 4.0 * pi / 3.0;
  const double k_factor = settings.vv10_b * 1.5 * pi * std::pow( 9.0 * pi, -1.0 / 6.0 );
  const double beta_value = beta(settings);

  std::vector<double> omega(npts, 0.0), kappa(npts, 0.0), rp_weight(npts, 0.0);
  std::vector<double> domega_drho(npts, 0.0), domega_dgamma(npts, 0.0), dkappa_drho(npts, 0.0);
  std::vector<char> active(npts, 0);

  for( std::size_t i = 0; i < npts; ++i ) {
    const double rho = grid.rho[i];
    if( rho < settings.vv10_tol ) continue;

    const double gamma = std::max( grid.gamma[i], 0.0 );
    const double rho2 = rho * rho;
    const double omega_tmp = settings.vv10_c * (gamma / rho2) * (gamma / rho2);
    const double omega_i = std::sqrt( omega_tmp + pi43 * rho );

    active[i] = 1;
    omega[i] = omega_i;
    kappa[i] = k_factor * std::pow( rho, 1.0 / 6.0 );
    rp_weight[i] = rho * grid.weights[i];
    domega_drho[i] = (0.5 * pi43 * rho - 2.0 * omega_tmp) / omega_i;
    domega_dgamma[i] = gamma > 0.0 ? omega_tmp * rho / (gamma * omega_i) : 0.0;
    dkappa_drho[i] = kappa[i] / 6.0;
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    double f = 0.0;
    double u = 0.0;
    double w = 0.0;
    const double* coord_i = grid.coords + 3*i;

    for( std::size_t j = 0; j < npts; ++j ) {
      if( not active[j] ) continue;

      const double* coord_j = grid.coords + 3*j;
      const double dx = coord_j[0] - coord_i[0];
      const double dy = coord_j[1] - coord_i[1];
      const double dz = coord_j[2] - coord_i[2];
      const double r2 = dx*dx + dy*dy + dz*dz;
      const double gp = r2 * omega[j] + kappa[j];
      const double g  = r2 * omega[i] + kappa[i];
      const double gt = g + gp;
      double t = rp_weight[j] / (g * gp * gt);

      f += t;
      t *= 1.0/g + 1.0/gt;
      u += t;
      w += t * r2;
    }

    f *= -1.5;
    corr.eps[i]    += beta_value + 0.5 * f;
    corr.vrho[i]   += beta_value + f + 1.5 * (u * dkappa_drho[i] + w * domega_drho[i]);
    corr.vgamma[i] += 1.5 * w * domega_dgamma[i];
  }
}

inline void eval_fxc_rks( const IntegratorSettingsNLC& settings,
                          ResponseGridView grid, TrialView trial,
                          ResponseCorrectionsView corr ) {
  if( grid.npts == 0 ) return;

  const auto npts = grid.npts;
  const double pi = 3.141592653589793238462643383279502884;
  const double pi43 = 4.0 * pi / 3.0;
  const double k_factor = settings.vv10_b * 1.5 * pi * std::pow( 9.0 * pi, -1.0 / 6.0 );

  std::vector<double> omega(npts, 0.0), kappa(npts, 0.0);
  std::vector<double> uvec(npts, 0.0), wvec(npts, 0.0), avec(npts, 0.0), bvec(npts, 0.0), cvec(npts, 0.0);
  std::vector<double> domega_drho(npts, 0.0), domega_dgamma(npts, 0.0), dkappa_drho(npts, 0.0);
  std::vector<double> d2omega_drho2(npts, 0.0), d2omega_dgamma2(npts, 0.0), d2omega_drho_dgamma(npts, 0.0);
  std::vector<double> d2kappa_drho2(npts, 0.0), gamma_t(npts, 0.0), f_rho_t(npts, 0.0), f_gamma_t(npts, 0.0);
  std::vector<char> active(npts, 0);

  for( std::size_t i = 0; i < npts; ++i ) {
    const double rho = grid.rho[i];
    if( rho < settings.vv10_tol ) continue;

    const double gamma = std::max( grid.gamma[i], 0.0 );
    const double rho_inv = 1.0 / rho;
    const double rho2_inv = rho_inv * rho_inv;
    const double rho3_inv = rho2_inv * rho_inv;
    const double rho4_inv = rho2_inv * rho2_inv;
    const double rho5_inv = rho4_inv * rho_inv;
    const double gamma2 = gamma * gamma;
    const double omega2 = settings.vv10_c * gamma2 * rho4_inv + pi43 * rho;
    const double omega_i = std::sqrt( omega2 );
    const double omega_inv = 1.0 / omega_i;
    const double omega3_inv = omega_inv / omega2;

    active[i] = 1;
    omega[i] = omega_i;
    kappa[i] = k_factor * std::pow( rho, 1.0 / 6.0 );
    domega_drho[i] = 0.5 * (pi43 - 4.0 * settings.vv10_c * gamma2 * rho5_inv) * omega_inv;
    domega_dgamma[i] = settings.vv10_c * gamma * rho4_inv * omega_inv;
    dkappa_drho[i] = k_factor * (1.0 / 6.0) * std::pow( rho, -5.0 / 6.0 );
    d2omega_drho2[i] = (-0.25 * pi43 * pi43
      + 12.0 * pi43 * settings.vv10_c * gamma2 * rho5_inv
      + 6.0 * settings.vv10_c * settings.vv10_c * gamma2 * gamma2 * rho5_inv * rho5_inv) * omega3_inv;
    d2omega_dgamma2[i] = pi43 * settings.vv10_c * rho3_inv * omega3_inv;
    d2omega_drho_dgamma[i] = -settings.vv10_c * gamma *
      (4.5 * pi43 * rho4_inv + 2.0 * settings.vv10_c * gamma2 * rho4_inv * rho5_inv) * omega3_inv;
    d2kappa_drho2[i] = k_factor * (-5.0 / 36.0) * std::pow( rho, -11.0 / 6.0 );
    gamma_t[i] = 2.0 * (grid.grad_x[i] * trial.grad_x[i]
      + grid.grad_y[i] * trial.grad_y[i] + grid.grad_z[i] * trial.grad_z[i]);
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    const double omega_i = omega[i];
    const double kappa_i = kappa[i];
    const double* coord_i = grid.coords + 3*i;

    double u_i = 0.0;
    double w_i = 0.0;
    double a_i = 0.0;
    double b_i = 0.0;
    double c_i = 0.0;

    for( std::size_t j = 0; j < npts; ++j ) {
      if( not active[j] ) continue;

      const double* coord_j = grid.coords + 3*j;
      const double dx = coord_i[0] - coord_j[0];
      const double dy = coord_i[1] - coord_j[1];
      const double dz = coord_i[2] - coord_j[2];
      const double r2 = dx*dx + dy*dy + dz*dz;
      const double g_ij = omega_i * r2 + kappa_i;
      const double g_ji = omega[j] * r2 + kappa[j];
      const double g_ij_inv = 1.0 / g_ij;
      const double g_sum_inv = 1.0 / (g_ij + g_ji);
      const double phi_ij = -1.5 / g_ji * g_ij_inv * g_sum_inv;
      const double e_ij = grid.weights[j] * grid.rho[j] * phi_ij;
      const double u_ij = e_ij * (g_sum_inv + g_ij_inv);
      const double w_ij = u_ij * r2;
      const double a_ij = e_ij * (g_sum_inv * g_sum_inv + g_sum_inv * g_ij_inv + g_ij_inv * g_ij_inv);
      const double b_ij = a_ij * r2;
      const double c_ij = b_ij * r2;

      u_i += u_ij;
      w_i += w_ij;
      a_i += a_ij;
      b_i += b_ij;
      c_i += c_ij;
    }

    uvec[i] = -u_i;
    wvec[i] = -w_i;
    avec[i] = 2.0 * a_i;
    bvec[i] = 2.0 * b_i;
    cvec[i] = 2.0 * c_i;
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    const double omega_i = omega[i];
    const double kappa_i = kappa[i];
    const double rho_i = grid.rho[i];
    const double domega_drho_i = domega_drho[i];
    const double domega_dgamma_i = domega_dgamma[i];
    const double dkappa_drho_i = dkappa_drho[i];
    const double* coord_i = grid.coords + 3*i;

    double f_rho_t_i = 0.0;
    double f_gamma_t_i = 0.0;

    for( std::size_t j = 0; j < npts; ++j ) {
      if( not active[j] ) continue;

      const double* coord_j = grid.coords + 3*j;
      const double dx = coord_i[0] - coord_j[0];
      const double dy = coord_i[1] - coord_j[1];
      const double dz = coord_i[2] - coord_j[2];
      const double r2 = dx*dx + dy*dy + dz*dz;
      const double g_ij = omega_i * r2 + kappa_i;
      const double g_ji = omega[j] * r2 + kappa[j];
      const double g_ij_inv = 1.0 / g_ij;
      const double g_ji_inv = 1.0 / g_ji;
      const double g_sum_inv = 1.0 / (g_ij + g_ji);
      const double phi_ij = -1.5 * g_ij_inv * g_ji_inv * g_sum_inv;
      const double rho_dgdrho_i = rho_i * (r2 * domega_drho_i + dkappa_drho_i);
      const double rho_dgdrho_j = grid.rho[j] * (r2 * domega_drho[j] + dkappa_drho[j]);
      const double d2phi_dgij_dgji_over_phi = 2.0 * (g_sum_inv * g_sum_inv + g_ij_inv * g_ji_inv);

      const double f_rho_rho_ij = phi_ij * (rho_dgdrho_i * rho_dgdrho_j * d2phi_dgij_dgji_over_phi
        - rho_dgdrho_i * (g_sum_inv + g_ij_inv)
        - rho_dgdrho_j * (g_sum_inv + g_ji_inv) + 1.0);
      const double f_gamma_rho_ij = rho_i * domega_dgamma_i * r2 * phi_ij *
        (rho_dgdrho_j * d2phi_dgij_dgji_over_phi - (g_sum_inv + g_ij_inv));
      const double f_rho_gamma_ij = grid.rho[j] * domega_dgamma[j] * r2 * phi_ij *
        (rho_dgdrho_i * d2phi_dgij_dgji_over_phi - (g_sum_inv + g_ji_inv));
      const double f_gamma_gamma_ij = rho_i * grid.rho[j] * domega_dgamma_i * domega_dgamma[j]
        * r2 * r2 * phi_ij * d2phi_dgij_dgji_over_phi;

      f_rho_t_i += grid.weights[j] * (f_rho_rho_ij * trial.rho[j] + f_rho_gamma_ij * gamma_t[j]);
      f_gamma_t_i += grid.weights[j] * (f_gamma_rho_ij * trial.rho[j] + f_gamma_gamma_ij * gamma_t[j]);
    }

    const double f_rho_rho_ii = 2.0 * domega_drho_i * wvec[i] + 2.0 * dkappa_drho_i * uvec[i]
      + rho_i * (d2omega_drho2[i] * wvec[i] + d2kappa_drho2[i] * uvec[i]
      + dkappa_drho_i * dkappa_drho_i * avec[i] + domega_drho_i * domega_drho_i * cvec[i]
      + 2.0 * domega_drho_i * dkappa_drho_i * bvec[i]);
    const double f_gamma_rho_ii = domega_dgamma_i * wvec[i] + rho_i *
      (d2omega_drho_dgamma[i] * wvec[i] + domega_dgamma_i *
      (dkappa_drho_i * bvec[i] + domega_drho_i * cvec[i]));
    const double f_gamma_gamma_ii = rho_i *
      (d2omega_dgamma2[i] * wvec[i] + domega_dgamma_i * domega_dgamma_i * cvec[i]);

    f_rho_t_i += f_rho_rho_ii * trial.rho[i] + f_gamma_rho_ii * gamma_t[i];
    f_gamma_t_i += f_gamma_rho_ii * trial.rho[i] + f_gamma_gamma_ii * gamma_t[i];

    f_rho_t[i] = f_rho_t_i;
    f_gamma_t[i] = f_gamma_t_i;
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    const double f_gamma = grid.rho[i] * domega_dgamma[i] * wvec[i];
    corr.A[i] += f_rho_t[i];
    corr.B[3*i]   += 2.0 * (grid.grad_x[i] * f_gamma_t[i] + trial.grad_x[i] * f_gamma);
    corr.B[3*i+1] += 2.0 * (grid.grad_y[i] * f_gamma_t[i] + trial.grad_y[i] * f_gamma);
    corr.B[3*i+2] += 2.0 * (grid.grad_z[i] * f_gamma_t[i] + trial.grad_z[i] * f_gamma);
  }
}

inline void eval_grid_gradient( const IntegratorSettingsNLC& settings,
                                GridView grid, GridGradientView grad ) {
  if( grid.npts == 0 ) return;

  const auto npts = grid.npts;
  const double pi = 3.141592653589793238462643383279502884;
  const double pi43 = 4.0 * pi / 3.0;
  const double k_factor = settings.vv10_b * 1.5 * pi * std::pow( 9.0 * pi, -1.0 / 6.0 );

  std::vector<double> omega(npts, 0.0), kappa(npts, 0.0), rho_weight(npts, 0.0);
  std::vector<char> active(npts, 0);

  for( std::size_t i = 0; i < npts; ++i ) {
    const double rho = grid.rho[i];
    if( rho < settings.vv10_tol ) continue;

    const double gamma = std::max( grid.gamma[i], 0.0 );
    const double rho2 = rho * rho;
    const double omega_tmp = settings.vv10_c * (gamma / rho2) * (gamma / rho2);
    active[i] = 1;
    omega[i] = std::sqrt( omega_tmp + pi43 * rho );
    kappa[i] = k_factor * std::pow( rho, 1.0 / 6.0 );
    rho_weight[i] = rho * grid.weights[i];
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    const double* coord_i = grid.coords + 3*i;

    for( std::size_t j = 0; j < npts; ++j ) {
      if( not active[j] ) continue;

      const double* coord_j = grid.coords + 3*j;
      const double dx = coord_j[0] - coord_i[0];
      const double dy = coord_j[1] - coord_i[1];
      const double dz = coord_j[2] - coord_i[2];
      const double r2 = dx*dx + dy*dy + dz*dz;
      const double gp = r2 * omega[j] + kappa[j];
      const double g  = r2 * omega[i] + kappa[i];
      const double gt = g + gp;
      const double t = rho_weight[j] / (g * gp * gt);
      const double q = t * (omega[i]/g + omega[j]/gp + (omega[i]+omega[j])/gt);

      fx += q * dx;
      fy += q * dy;
      fz += q * dz;
    }

    grad.x[i] += -3.0 * fx;
    grad.y[i] += -3.0 * fy;
    grad.z[i] += -3.0 * fz;
  }
}

inline void eval_grid_gradient_excluding_same_parent( const IntegratorSettingsNLC& settings,
                                                      GridView grid,
                                                      const std::vector<unsigned long long>& parent,
                                                      GridGradientView grad ) {
  if( grid.npts == 0 ) return;

  const auto npts = grid.npts;
  const double pi = 3.141592653589793238462643383279502884;
  const double pi43 = 4.0 * pi / 3.0;
  const double k_factor = settings.vv10_b * 1.5 * pi * std::pow( 9.0 * pi, -1.0 / 6.0 );

  std::vector<double> omega(npts, 0.0), kappa(npts, 0.0), rho_weight(npts, 0.0);
  std::vector<char> active(npts, 0);

  for( std::size_t i = 0; i < npts; ++i ) {
    const double rho = grid.rho[i];
    if( rho < settings.vv10_tol ) continue;

    const double gamma = std::max( grid.gamma[i], 0.0 );
    const double rho2 = rho * rho;
    const double omega_tmp = settings.vv10_c * (gamma / rho2) * (gamma / rho2);
    active[i] = 1;
    omega[i] = std::sqrt( omega_tmp + pi43 * rho );
    kappa[i] = k_factor * std::pow( rho, 1.0 / 6.0 );
    rho_weight[i] = rho * grid.weights[i];
  }

  for( std::size_t i = 0; i < npts; ++i ) {
    if( not active[i] ) continue;

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    const double* coord_i = grid.coords + 3*i;

    for( std::size_t j = 0; j < npts; ++j ) {
      if( not active[j] or parent[i] == parent[j] ) continue;

      const double* coord_j = grid.coords + 3*j;
      const double dx = coord_j[0] - coord_i[0];
      const double dy = coord_j[1] - coord_i[1];
      const double dz = coord_j[2] - coord_i[2];
      const double r2 = dx*dx + dy*dy + dz*dz;
      const double gp = r2 * omega[j] + kappa[j];
      const double g  = r2 * omega[i] + kappa[i];
      const double gt = g + gp;
      const double t = rho_weight[j] / (g * gp * gt);
      const double q = t * (omega[i]/g + omega[j]/gp + (omega[i]+omega[j])/gt);

      fx += q * dx;
      fy += q * dy;
      fz += q * dz;
    }

    grad.x[i] += -3.0 * fx;
    grad.y[i] += -3.0 * fy;
    grad.z[i] += -3.0 * fz;
  }
}

} // namespace GauXC::detail::vv10