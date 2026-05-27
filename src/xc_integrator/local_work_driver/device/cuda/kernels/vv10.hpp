#pragma once

#include <gauxc/gauxc_config.hpp>

#ifdef GAUXC_HAS_CUDA

#include "integrator_util/vv10.hpp"
#include <cuda_runtime.h>
#include <cstdint>

namespace GauXC::detail::vv10 {

void eval_omega_kappa_cuda(cudaStream_t stream, std::int32_t npts,
  const double* rho, const double* gamma, Parameters params, double* omega,
  double* kappa, double* domega_drho, double* domega_dgamma,
  double* d2omega_drho2, double* d2omega_dgamma2,
  double* d2omega_drho_dgamma, double* dkappa_drho,
  double* d2kappa_drho2);

void eval_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* omega, const double* kappa,
  std::int32_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* F, double* U, double* W);

void eval_exc_vxc_from_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* rho, const double* gamma, const double* F, const double* U,
  const double* W, Parameters params, double* eps, double* vrho,
  double* vgamma);

void eval_exc_vxc_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* rho, const double* gamma,
  const double* omega, const double* kappa, const double* rho_weight,
  Parameters params, double* eps, double* vrho, double* vgamma);

void eval_grid_gradient_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* omega, const double* kappa,
  std::int32_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* gradient);

void eval_grid_gradient_excluding_same_parent_cuda(cudaStream_t stream,
  std::int32_t npts, const double* coords, const double* omega,
  const double* kappa, const double* rho_weight,
  const std::int32_t* parent, double* gradient);

void eval_hessian_intermediates_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa, double* U, double* W,
  double* A, double* B, double* C, double* E);

void eval_hessian_contraction_cuda(cudaStream_t stream, std::int32_t npts,
  std::int32_t ntrial, const double* coords, const double* weights,
  const double* rho, const double* omega, const double* kappa,
  const double* U, const double* W, const double* A, const double* B,
  const double* C, const double* domega_drho,
  const double* domega_dgamma, const double* dkappa_drho,
  const double* d2omega_drho2, const double* d2omega_dgamma2,
  const double* d2omega_drho_dgamma, const double* d2kappa_drho2,
  const double* rho_t, const double* gamma_t, double* f_rho_t,
  double* f_gamma_t);

void eval_grid_response_cuda(cudaStream_t stream, std::int32_t npts,
  std::int32_t natoms, const double* coords, const double* weights,
  const double* rho, const double* omega, const double* kappa,
  const std::int32_t* grid_associated_atom, double* Egr, double* Ugr,
  double* Wgr);

}

#endif