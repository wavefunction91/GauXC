#include "integrator_util/vv10.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace GauXC::detail::vv10 {

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double four_pi_over_three = 4.0 * pi / 3.0;

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

Vec3 load_point(const double* coords, std::int64_t i) {
  return {coords[3*i + 0], coords[3*i + 1], coords[3*i + 2]};
}

double distance2(Vec3 a, Vec3 b) {
  const auto dx = a.x - b.x;
  const auto dy = a.y - b.y;
  const auto dz = a.z - b.z;
  return dx*dx + dy*dy + dz*dz;
}

}

double kappa_prefactor(const Parameters& params) {
  return params.b * 1.5 * pi * std::pow(9.0 * pi, -1.0/6.0);
}

double beta(const Parameters& params) {
  return std::pow(3.0 / (params.b * params.b), 0.75) / 32.0;
}

void eval_omega_kappa(std::int64_t npts, const double* rho, const double* gamma,
  const Parameters& params, double* omega, double* kappa,
  double* domega_drho, double* domega_dgamma, double* d2omega_drho2,
  double* d2omega_dgamma2, double* d2omega_drho_dgamma,
  double* dkappa_drho, double* d2kappa_drho2) {

  const auto kappa_pref = kappa_prefactor(params);

  #pragma omp parallel for schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    const auto rho_i = rho[i];
    const auto gamma_i = gamma[i];

    if(rho_i < params.rho_threshold) {
      omega[i] = 0.0;
      kappa[i] = 0.0;
      if(domega_drho)          domega_drho[i] = 0.0;
      if(domega_dgamma)        domega_dgamma[i] = 0.0;
      if(d2omega_drho2)        d2omega_drho2[i] = 0.0;
      if(d2omega_dgamma2)      d2omega_dgamma2[i] = 0.0;
      if(d2omega_drho_dgamma) d2omega_drho_dgamma[i] = 0.0;
      if(dkappa_drho)          dkappa_drho[i] = 0.0;
      if(d2kappa_drho2)        d2kappa_drho2[i] = 0.0;
      continue;
    }

    const auto rho_1 = 1.0 / rho_i;
    const auto rho_2 = rho_1 * rho_1;
    const auto rho_3 = rho_1 * rho_2;
    const auto rho_4 = rho_2 * rho_2;
    const auto rho_5 = rho_1 * rho_4;
    const auto gamma2 = gamma_i * gamma_i;
    const auto omega2 = params.c * gamma2 * rho_4 + four_pi_over_three * rho_i;
    const auto omega_i = std::sqrt(omega2);
    const auto omega_1 = 1.0 / omega_i;

    omega[i] = omega_i;
    kappa[i] = kappa_pref * std::pow(rho_i, 1.0/6.0);

    if(domega_drho) {
      domega_drho[i] = 0.5 * (four_pi_over_three - 4.0 * params.c * gamma2 * rho_5) * omega_1;
    }
    if(domega_dgamma) {
      domega_dgamma[i] = params.c * gamma_i * rho_4 * omega_1;
    }
    if(d2omega_drho2 or d2omega_dgamma2 or d2omega_drho_dgamma) {
      const auto omega_3 = omega_1 / omega2;
      if(d2omega_drho2) {
        d2omega_drho2[i] = (-0.25 * four_pi_over_three * four_pi_over_three
          + 12.0 * four_pi_over_three * params.c * gamma2 * rho_5
          + 6.0 * params.c * params.c * gamma2 * gamma2 * rho_5 * rho_5) * omega_3;
      }
      if(d2omega_dgamma2) {
        d2omega_dgamma2[i] = four_pi_over_three * params.c * rho_3 * omega_3;
      }
      if(d2omega_drho_dgamma) {
        d2omega_drho_dgamma[i] = -params.c * gamma_i * (4.5 * four_pi_over_three * rho_4
          + 2.0 * params.c * gamma2 * rho_4 * rho_5) * omega_3;
      }
    }
    if(dkappa_drho) {
      dkappa_drho[i] = kappa_pref * (1.0/6.0) * std::pow(rho_i, -5.0/6.0);
    }
    if(d2kappa_drho2) {
      d2kappa_drho2[i] = kappa_pref * (-5.0/36.0) * std::pow(rho_i, -11.0/6.0);
    }
  }
}

void eval_fock(std::int64_t npts, const double* coords, const double* omega,
  const double* kappa, std::int64_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* F, double* U, double* W) {

  #pragma omp parallel for schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    const auto r_i = load_point(coords, i);
    const auto omega_i = omega[i];
    const auto kappa_i = kappa[i];
    double F_i = 0.0;
    double U_i = 0.0;
    double W_i = 0.0;

    for(std::int64_t j = 0; j < vv_npts; ++j) {
      const auto r_j = load_point(vv_coords, j);
      const auto r_ij2 = distance2(r_j, r_i);
      const auto g_ij = r_ij2 * omega_i + kappa_i;
      const auto g_ji = r_ij2 * vv_omega[j] + vv_kappa[j];
      const auto g_sum = g_ij + g_ji;
      auto pair = vv_rho_weight[j] / (g_ij * g_ji * g_sum);
      F_i += pair;
      pair *= 1.0/g_ij + 1.0/g_sum;
      U_i += pair;
      W_i += pair * r_ij2;
    }

    F[i] = -1.5 * F_i;
    U[i] = U_i;
    W[i] = W_i;
  }
}

void eval_exc_vxc(std::int64_t npts, const double* coords,
  const double* weights, const double* rho, const double* gamma,
  const Parameters& params, double* eps, double* vrho, double* vgamma) {

  std::vector<double> omega(npts), kappa(npts), domega_drho(npts), domega_dgamma(npts);
  std::vector<double> F(npts), U(npts), W(npts), rho_weight(npts);

  eval_omega_kappa(npts, rho, gamma, params, omega.data(), kappa.data(),
    domega_drho.data(), domega_dgamma.data());

  for(std::int64_t i = 0; i < npts; ++i) {
    rho_weight[i] = rho[i] >= params.rho_threshold ? rho[i] * weights[i] : 0.0;
  }

  eval_fock(npts, coords, omega.data(), kappa.data(), npts, coords,
    rho_weight.data(), omega.data(), kappa.data(), F.data(), U.data(), W.data());

  const auto beta_value = beta(params);
  const auto kappa_pref = kappa_prefactor(params);

  #pragma omp parallel for schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    if(rho[i] < params.rho_threshold) {
      eps[i] = 0.0;
      vrho[i] = 0.0;
      vgamma[i] = 0.0;
      continue;
    }

    const auto scaled_dkappa_drho = (1.0/6.0) * kappa_pref * std::pow(rho[i], 1.0/6.0);
    const auto scaled_domega_drho = rho[i] * domega_drho[i];
    const auto scaled_domega_dgamma = rho[i] * domega_dgamma[i];
    eps[i] = beta_value + 0.5 * F[i];
    vrho[i] = beta_value + F[i] + 1.5 * (U[i] * scaled_dkappa_drho + W[i] * scaled_domega_drho);
    vgamma[i] = 1.5 * W[i] * scaled_domega_dgamma;
  }
}

void eval_grid_gradient(std::int64_t npts, const double* coords,
  const double* omega, const double* kappa, std::int64_t vv_npts,
  const double* vv_coords, const double* vv_rho_weight,
  const double* vv_omega, const double* vv_kappa, double* gradient) {

  #pragma omp parallel for schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    const auto r_i = load_point(coords, i);
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;

    for(std::int64_t j = 0; j < vv_npts; ++j) {
      const auto r_j = load_point(vv_coords, j);
      const auto dx = r_j.x - r_i.x;
      const auto dy = r_j.y - r_i.y;
      const auto dz = r_j.z - r_i.z;
      const auto r_ij2 = dx*dx + dy*dy + dz*dz;
      const auto g_ij = r_ij2 * omega[i] + kappa[i];
      const auto g_ji = r_ij2 * vv_omega[j] + vv_kappa[j];
      const auto g_sum = g_ij + g_ji;
      const auto pair = vv_rho_weight[j] / (g_ij * g_ji * g_sum);
      const auto prefactor = pair * (omega[i]/g_ij + vv_omega[j]/g_ji + (omega[i] + vv_omega[j])/g_sum);
      gx += prefactor * dx;
      gy += prefactor * dy;
      gz += prefactor * dz;
    }

    gradient[3*i + 0] = -3.0 * gx;
    gradient[3*i + 1] = -3.0 * gy;
    gradient[3*i + 2] = -3.0 * gz;
  }
}

void eval_hessian_intermediates(std::int64_t npts, const double* coords,
  const double* weights, const double* rho, const double* omega,
  const double* kappa, double* U, double* W, double* A, double* B,
  double* C, double* E) {

  #pragma omp parallel for schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    const auto r_i = load_point(coords, i);
    double U_i = 0.0;
    double W_i = 0.0;
    double A_i = 0.0;
    double B_i = 0.0;
    double C_i = 0.0;
    double E_i = 0.0;

    for(std::int64_t j = 0; j < npts; ++j) {
      const auto r_j = load_point(coords, j);
      const auto r_ij2 = distance2(r_i, r_j);
      const auto g_ij = omega[i] * r_ij2 + kappa[i];
      const auto g_ji = omega[j] * r_ij2 + kappa[j];
      const auto g_ij_1 = 1.0 / g_ij;
      const auto g_ji_1 = 1.0 / g_ji;
      const auto g_sum_1 = 1.0 / (g_ij + g_ji);
      const auto phi = -1.5 * g_ij_1 * g_ji_1 * g_sum_1;
      const auto E_ij = weights[j] * rho[j] * phi;
      const auto U_ij = E_ij * (g_sum_1 + g_ij_1);
      const auto W_ij = U_ij * r_ij2;
      const auto A_ij = E_ij * (g_sum_1*g_sum_1 + g_sum_1*g_ij_1 + g_ij_1*g_ij_1);
      const auto B_ij = A_ij * r_ij2;
      const auto C_ij = B_ij * r_ij2;
      U_i += U_ij;
      W_i += W_ij;
      A_i += A_ij;
      B_i += B_ij;
      C_i += C_ij;
      E_i += E_ij;
    }

    U[i] = -U_i;
    W[i] = -W_i;
    A[i] = 2.0 * A_i;
    B[i] = 2.0 * B_i;
    C[i] = 2.0 * C_i;
    E[i] = E_i;
  }
}

void eval_hessian_contraction(std::int64_t npts, std::int64_t ntrial,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa, const double* U,
  const double* W, const double* A, const double* B, const double* C,
  const double* domega_drho, const double* domega_dgamma,
  const double* dkappa_drho, const double* d2omega_drho2,
  const double* d2omega_dgamma2, const double* d2omega_drho_dgamma,
  const double* d2kappa_drho2, const double* rho_t,
  const double* gamma_t, double* f_rho_t, double* f_gamma_t) {

  #pragma omp parallel for collapse(2) schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    for(std::int64_t trial = 0; trial < ntrial; ++trial) {
      const auto r_i = load_point(coords, i);
      double f_rho_i = 0.0;
      double f_gamma_i = 0.0;

      for(std::int64_t j = 0; j < npts; ++j) {
        const auto r_j = load_point(coords, j);
        const auto r_ij2 = distance2(r_i, r_j);
        const auto g_ij = omega[i] * r_ij2 + kappa[i];
        const auto g_ji = omega[j] * r_ij2 + kappa[j];
        const auto g_ij_1 = 1.0 / g_ij;
        const auto g_ji_1 = 1.0 / g_ji;
        const auto g_sum_1 = 1.0 / (g_ij + g_ji);
        const auto phi = -1.5 * g_ij_1 * g_ji_1 * g_sum_1;
        const auto rho_dgdrho_i = rho[i] * (r_ij2 * domega_drho[i] + dkappa_drho[i]);
        const auto rho_dgdrho_j = rho[j] * (r_ij2 * domega_drho[j] + dkappa_drho[j]);
        const auto d2phi_dgij_dgji_over_phi = 2.0 * (g_sum_1*g_sum_1 + g_ij_1*g_ji_1);

        const auto f_rho_rho_ij = phi * (rho_dgdrho_i * rho_dgdrho_j * d2phi_dgij_dgji_over_phi
          - rho_dgdrho_i * (g_sum_1 + g_ij_1) - rho_dgdrho_j * (g_sum_1 + g_ji_1) + 1.0);
        const auto f_gamma_rho_ij = rho[i] * domega_dgamma[i] * r_ij2 * phi
          * (rho_dgdrho_j * d2phi_dgij_dgji_over_phi - (g_sum_1 + g_ij_1));
        const auto f_rho_gamma_ij = rho[j] * domega_dgamma[j] * r_ij2 * phi
          * (rho_dgdrho_i * d2phi_dgij_dgji_over_phi - (g_sum_1 + g_ji_1));
        const auto f_gamma_gamma_ij = rho[i] * rho[j] * domega_dgamma[i] * domega_dgamma[j]
          * r_ij2 * r_ij2 * phi * d2phi_dgij_dgji_over_phi;

        const auto offset = trial * npts + j;
        f_rho_i += weights[j] * (f_rho_rho_ij * rho_t[offset] + f_rho_gamma_ij * gamma_t[offset]);
        f_gamma_i += weights[j] * (f_gamma_rho_ij * rho_t[offset] + f_gamma_gamma_ij * gamma_t[offset]);
      }

      const auto f_rho_rho_ii = 2.0 * domega_drho[i] * W[i] + 2.0 * dkappa_drho[i] * U[i]
        + rho[i] * (d2omega_drho2[i] * W[i] + d2kappa_drho2[i] * U[i]
        + dkappa_drho[i] * dkappa_drho[i] * A[i] + domega_drho[i] * domega_drho[i] * C[i]
        + 2.0 * domega_drho[i] * dkappa_drho[i] * B[i]);
      const auto f_gamma_rho_ii = domega_dgamma[i] * W[i] + rho[i] * (d2omega_drho_dgamma[i] * W[i]
        + domega_dgamma[i] * (dkappa_drho[i] * B[i] + domega_drho[i] * C[i]));
      const auto f_gamma_gamma_ii = rho[i] * (d2omega_dgamma2[i] * W[i]
        + domega_dgamma[i] * domega_dgamma[i] * C[i]);
      const auto ioffset = trial * npts + i;

      f_rho_i += f_rho_rho_ii * rho_t[ioffset] + f_gamma_rho_ii * gamma_t[ioffset];
      f_gamma_i += f_gamma_rho_ii * rho_t[ioffset] + f_gamma_gamma_ii * gamma_t[ioffset];

      f_rho_t[ioffset] = f_rho_i;
      f_gamma_t[ioffset] = f_gamma_i;
    }
  }
}

void eval_grid_response(std::int64_t npts, std::int64_t natoms,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa,
  const std::int32_t* grid_associated_atom, double* Egr, double* Ugr,
  double* Wgr) {

  #pragma omp parallel for collapse(2) schedule(static)
  for(std::int64_t i = 0; i < npts; ++i) {
    for(std::int64_t atom = 0; atom < natoms; ++atom) {
      const auto i_atom = grid_associated_atom[i];
      const auto base = atom * 3 * npts + i;
      if(i_atom < 0) {
        Egr[base + 0*npts] = 0.0;
        Egr[base + 1*npts] = 0.0;
        Egr[base + 2*npts] = 0.0;
        Ugr[base + 0*npts] = 0.0;
        Ugr[base + 1*npts] = 0.0;
        Ugr[base + 2*npts] = 0.0;
        Wgr[base + 0*npts] = 0.0;
        Wgr[base + 1*npts] = 0.0;
        Wgr[base + 2*npts] = 0.0;
        continue;
      }

      const auto r_i = load_point(coords, i);
      const auto i_in_atom = i_atom == atom;
      Vec3 Eg{};
      Vec3 Ug{};
      Vec3 Wg{};

      for(std::int64_t j = 0; j < npts; ++j) {
        const auto j_atom = grid_associated_atom[j];
        if(j_atom < 0) continue;
        const auto j_in_atom = j_atom == atom;
        if((not i_in_atom and not j_in_atom) or (i_in_atom and j_in_atom)) continue;

        const auto r_j = load_point(coords, j);
        const Vec3 r_ji{r_j.x - r_i.x, r_j.y - r_i.y, r_j.z - r_i.z};
        const auto r_ij2 = r_ji.x*r_ji.x + r_ji.y*r_ji.y + r_ji.z*r_ji.z;
        const auto g_ij = omega[i] * r_ij2 + kappa[i];
        const auto g_ji = omega[j] * r_ij2 + kappa[j];
        const auto g_ij_1 = 1.0 / g_ij;
        const auto g_ji_1 = 1.0 / g_ji;
        const auto g_sum_1 = 1.0 / (g_ij + g_ji);
        const auto phi = -1.5 * g_ij_1 * g_ji_1 * g_sum_1;
        const auto E_ij = weights[j] * rho[j] * phi;
        const auto dphi_drj_over_phi = omega[i]*g_ij_1 + omega[j]*g_ji_1 + (omega[i] + omega[j])*g_sum_1;
        const auto d2phi_dgij_drj_over_phi = omega[i]*g_ij_1*g_ij_1 + (omega[i] + omega[j])*g_sum_1*g_sum_1;
        const auto dphi_dgij_over_phi = g_sum_1 + g_ij_1;
        const auto Egr_ij = E_ij * dphi_drj_over_phi;
        const auto Ugr_ij = E_ij * (dphi_drj_over_phi*dphi_dgij_over_phi + d2phi_dgij_drj_over_phi);
        const auto Wgr_ij = E_ij * (r_ij2 * (dphi_drj_over_phi*dphi_dgij_over_phi + d2phi_dgij_drj_over_phi) - dphi_dgij_over_phi);

        Eg.x += Egr_ij * r_ji.x;
        Eg.y += Egr_ij * r_ji.y;
        Eg.z += Egr_ij * r_ji.z;
        Ug.x += Ugr_ij * r_ji.x;
        Ug.y += Ugr_ij * r_ji.y;
        Ug.z += Ugr_ij * r_ji.z;
        Wg.x += Wgr_ij * r_ji.x;
        Wg.y += Wgr_ij * r_ji.y;
        Wg.z += Wgr_ij * r_ji.z;
      }

      const auto atom_sign = i_in_atom ? -1.0 : 1.0;
      Egr[base + 0*npts] = -2.0 * atom_sign * Eg.x;
      Egr[base + 1*npts] = -2.0 * atom_sign * Eg.y;
      Egr[base + 2*npts] = -2.0 * atom_sign * Eg.z;
      Ugr[base + 0*npts] =  2.0 * atom_sign * Ug.x;
      Ugr[base + 1*npts] =  2.0 * atom_sign * Ug.y;
      Ugr[base + 2*npts] =  2.0 * atom_sign * Ug.z;
      Wgr[base + 0*npts] =  2.0 * atom_sign * Wg.x;
      Wgr[base + 1*npts] =  2.0 * atom_sign * Wg.y;
      Wgr[base + 2*npts] =  2.0 * atom_sign * Wg.z;
    }
  }
}

}