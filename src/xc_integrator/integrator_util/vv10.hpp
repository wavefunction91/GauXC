#pragma once

#include <cstdint>

namespace GauXC::detail::vv10 {

struct Parameters {
  double b = 6.3;
  double c = 0.0093;
  double rho_threshold = 1e-8;
  std::int32_t pair_math_mode = 0;
};

double kappa_prefactor(const Parameters& params);
double beta(const Parameters& params);

void eval_omega_kappa(std::int64_t npts, const double* rho, const double* gamma,
  const Parameters& params, double* omega, double* kappa,
  double* domega_drho = nullptr, double* domega_dgamma = nullptr,
  double* d2omega_drho2 = nullptr, double* d2omega_dgamma2 = nullptr,
  double* d2omega_drho_dgamma = nullptr, double* dkappa_drho = nullptr,
  double* d2kappa_drho2 = nullptr);

void eval_fock(std::int64_t npts, const double* coords, const double* omega,
  const double* kappa, std::int64_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* F, double* U, double* W);

void eval_exc_vxc(std::int64_t npts, const double* coords,
  const double* weights, const double* rho, const double* gamma,
  const Parameters& params, double* eps, double* vrho, double* vgamma);

void eval_grid_gradient(std::int64_t npts, const double* coords,
  const double* omega, const double* kappa, std::int64_t vv_npts,
  const double* vv_coords, const double* vv_rho_weight,
  const double* vv_omega, const double* vv_kappa, double* gradient);

void eval_hessian_intermediates(std::int64_t npts, const double* coords,
  const double* weights, const double* rho, const double* omega,
  const double* kappa, double* U, double* W, double* A, double* B,
  double* C, double* E);

void eval_hessian_contraction(std::int64_t npts, std::int64_t ntrial,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa, const double* U,
  const double* W, const double* A, const double* B, const double* C,
  const double* domega_drho, const double* domega_dgamma,
  const double* dkappa_drho, const double* d2omega_drho2,
  const double* d2omega_dgamma2, const double* d2omega_drho_dgamma,
  const double* d2kappa_drho2, const double* rho_t,
  const double* gamma_t, double* f_rho_t, double* f_gamma_t);

void eval_grid_response(std::int64_t npts, std::int64_t natoms,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa,
  const std::int32_t* grid_associated_atom, double* Egr, double* Ugr,
  double* Wgr);

}