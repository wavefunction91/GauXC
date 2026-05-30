#include "vv10.hpp"

#ifdef GAUXC_HAS_CUDA

#include "exceptions/cuda_exception.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <cmath>

namespace GauXC::detail::vv10 {

namespace {

constexpr int threads_per_block = 128;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double four_pi_over_three = 4.0 * pi / 3.0;

__device__ double3 load_point(const double* coords, int i);
__device__ double distance2(double3 a, double3 b);

__device__ double cuda_kappa_prefactor(Parameters params) {
  return params.b * 1.5 * pi * pow(9.0 * pi, -1.0/6.0);
}

__device__ double cuda_beta(Parameters params) {
  return pow(3.0 / (params.b * params.b), 0.75) / 32.0;
}

__device__ double3 load_point(const double* coords, int i) {
  return {coords[3*i + 0], coords[3*i + 1], coords[3*i + 2]};
}

__device__ double distance2(double3 a, double3 b) {
  const auto dx = a.x - b.x;
  const auto dy = a.y - b.y;
  const auto dz = a.z - b.z;
  return dx*dx + dy*dy + dz*dz;
}

__global__ void omega_kappa_kernel(int npts, const double* rho,
  const double* gamma, Parameters params, double* omega, double* kappa,
  double* domega_drho, double* domega_dgamma, double* d2omega_drho2,
  double* d2omega_dgamma2, double* d2omega_drho_dgamma,
  double* dkappa_drho, double* d2kappa_drho2) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= npts) return;

  const auto rho_i = rho[i];
  const auto gamma_i = gamma[i];
  if(rho_i < params.rho_threshold) {
    omega[i] = 0.0;
    kappa[i] = 0.0;
    if(domega_drho) domega_drho[i] = 0.0;
    if(domega_dgamma) domega_dgamma[i] = 0.0;
    if(d2omega_drho2) d2omega_drho2[i] = 0.0;
    if(d2omega_dgamma2) d2omega_dgamma2[i] = 0.0;
    if(d2omega_drho_dgamma) d2omega_drho_dgamma[i] = 0.0;
    if(dkappa_drho) dkappa_drho[i] = 0.0;
    if(d2kappa_drho2) d2kappa_drho2[i] = 0.0;
    return;
  }

  const auto kappa_pref = cuda_kappa_prefactor(params);
  const auto rho_1 = 1.0 / rho_i;
  const auto rho_2 = rho_1 * rho_1;
  const auto rho_3 = rho_1 * rho_2;
  const auto rho_4 = rho_2 * rho_2;
  const auto rho_5 = rho_1 * rho_4;
  const auto gamma2 = gamma_i * gamma_i;
  const auto omega2 = params.c * gamma2 * rho_4 + four_pi_over_three * rho_i;
  const auto omega_i = sqrt(omega2);
  const auto omega_1 = 1.0 / omega_i;

  omega[i] = omega_i;
  kappa[i] = kappa_pref * pow(rho_i, 1.0/6.0);
  if(domega_drho) domega_drho[i] = 0.5 * (four_pi_over_three - 4.0 * params.c * gamma2 * rho_5) * omega_1;
  if(domega_dgamma) domega_dgamma[i] = params.c * gamma_i * rho_4 * omega_1;
  if(d2omega_drho2 || d2omega_dgamma2 || d2omega_drho_dgamma) {
    const auto omega_3 = omega_1 / omega2;
    if(d2omega_drho2) {
      d2omega_drho2[i] = (-0.25 * four_pi_over_three * four_pi_over_three
        + 12.0 * four_pi_over_three * params.c * gamma2 * rho_5
        + 6.0 * params.c * params.c * gamma2 * gamma2 * rho_5 * rho_5) * omega_3;
    }
    if(d2omega_dgamma2) d2omega_dgamma2[i] = four_pi_over_three * params.c * rho_3 * omega_3;
    if(d2omega_drho_dgamma) {
      d2omega_drho_dgamma[i] = -params.c * gamma_i * (4.5 * four_pi_over_three * rho_4
        + 2.0 * params.c * gamma2 * rho_4 * rho_5) * omega_3;
    }
  }
  if(dkappa_drho) dkappa_drho[i] = kappa_pref * (1.0/6.0) * pow(rho_i, -5.0/6.0);
  if(d2kappa_drho2) d2kappa_drho2[i] = kappa_pref * (-5.0/36.0) * pow(rho_i, -11.0/6.0);
}

__global__ void fock_kernel(int npts, const double* coords,
  const double* omega, const double* kappa, int vv_npts,
  const double* vv_coords, const double* vv_rho_weight,
  const double* vv_omega, const double* vv_kappa, double* F, double* U,
  double* W) {

  __shared__ double3 coords_tile[threads_per_block];
  __shared__ double rho_weight_tile[threads_per_block];
  __shared__ double omega_tile[threads_per_block];
  __shared__ double kappa_tile[threads_per_block];

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  const auto active_i = i < npts;

  const auto r_i = active_i ? load_point(coords, i) : double3{0.0, 0.0, 0.0};
  const auto omega_i = active_i ? omega[i] : 0.0;
  const auto kappa_i = active_i ? kappa[i] : 0.0;
  double F_i = 0.0;
  double U_i = 0.0;
  double W_i = 0.0;

  for(int tile = 0; tile < vv_npts; tile += blockDim.x) {
    const auto j_load = tile + threadIdx.x;
    if(j_load < vv_npts) {
      coords_tile[threadIdx.x] = load_point(vv_coords, j_load);
      rho_weight_tile[threadIdx.x] = vv_rho_weight[j_load];
      omega_tile[threadIdx.x] = vv_omega[j_load];
      kappa_tile[threadIdx.x] = vv_kappa[j_load];
    }
    __syncthreads();

    if(active_i) {
      const auto tile_npts = min(blockDim.x, vv_npts - tile);
      for(int j = 0; j < tile_npts; ++j) {
        const auto r_ij2 = distance2(coords_tile[j], r_i);
        const auto g_ij = r_ij2 * omega_i + kappa_i;
        const auto g_ji = r_ij2 * omega_tile[j] + kappa_tile[j];
        const auto g_sum = g_ij + g_ji;
        const auto g_ij_inv = 1.0 / g_ij;
        const auto g_ji_inv = 1.0 / g_ji;
        const auto g_sum_inv = 1.0 / g_sum;
        auto pair = rho_weight_tile[j] * g_ij_inv * g_ji_inv * g_sum_inv;
        F_i += pair;
        pair *= g_ij_inv + g_sum_inv;
        U_i += pair;
        W_i += pair * r_ij2;
      }
    }
    __syncthreads();
  }

  if(active_i) {
    F[i] = -1.5 * F_i;
    U[i] = U_i;
    W[i] = W_i;
  }
}

__global__ void exc_vxc_from_fock_kernel(int npts, const double* rho,
  const double* gamma, const double* F, const double* U, const double* W,
  Parameters params, double* eps, double* vrho, double* vgamma) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= npts) return;
  if(rho[i] < params.rho_threshold) {
    eps[i] = 0.0;
    vrho[i] = 0.0;
    vgamma[i] = 0.0;
    return;
  }

  const auto rho_i = rho[i];
  const auto gamma_i = gamma[i];
  const auto rho_1 = 1.0 / rho_i;
  const auto rho_2 = rho_1 * rho_1;
  const auto rho_4 = rho_2 * rho_2;
  const auto rho_5 = rho_1 * rho_4;
  const auto gamma2 = gamma_i * gamma_i;
  const auto omega2 = params.c * gamma2 * rho_4 + four_pi_over_three * rho_i;
  const auto omega_i = sqrt(omega2);
  const auto domega_drho = 0.5 * (four_pi_over_three - 4.0 * params.c * gamma2 * rho_5) / omega_i;
  const auto domega_dgamma = params.c * gamma_i * rho_4 / omega_i;
  const auto kappa_pref = cuda_kappa_prefactor(params);
  const auto scaled_dkappa_drho = (1.0/6.0) * kappa_pref * pow(rho_i, 1.0/6.0);

  eps[i] = cuda_beta(params) + 0.5 * F[i];
  vrho[i] = cuda_beta(params) + F[i] + 1.5 * (U[i] * scaled_dkappa_drho + W[i] * rho_i * domega_drho);
  vgamma[i] = 1.5 * W[i] * rho_i * domega_dgamma;
}

template <bool UseFloatPair>
__global__ void exc_vxc_fock_kernel(int npts, const double* coords,
  const double* rho, const double* gamma, const double* omega,
  const double* kappa, const double* rho_weight, Parameters params,
  double* eps, double* vrho, double* vgamma) {

  __shared__ double3 coords_tile[threads_per_block];
  __shared__ double rho_weight_tile[threads_per_block];
  __shared__ double omega_tile[threads_per_block];
  __shared__ double kappa_tile[threads_per_block];
  __shared__ float rho_weight_tile_f[threads_per_block];
  __shared__ float omega_tile_f[threads_per_block];
  __shared__ float kappa_tile_f[threads_per_block];

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  const auto active_i = i < npts;

  const auto rho_i = active_i ? rho[i] : 0.0;
  const auto active_rho_i = active_i && rho_i >= params.rho_threshold;
  const auto r_i = active_rho_i ? load_point(coords, i) : double3{0.0, 0.0, 0.0};
  const auto omega_i = active_rho_i ? omega[i] : 0.0;
  const auto kappa_i = active_rho_i ? kappa[i] : 0.0;

  double F_i = 0.0;
  double U_i = 0.0;
  double W_i = 0.0;

  const float omega_i_f = static_cast<float>(omega_i);
  const float kappa_i_f = static_cast<float>(kappa_i);

  for(int tile = 0; tile < npts; tile += blockDim.x) {
    const auto j_load = tile + threadIdx.x;
    if(j_load < npts) {
      coords_tile[threadIdx.x] = load_point(coords, j_load);
      rho_weight_tile[threadIdx.x] = rho_weight[j_load];
      omega_tile[threadIdx.x] = omega[j_load];
      kappa_tile[threadIdx.x] = kappa[j_load];
      if constexpr(UseFloatPair) {
        rho_weight_tile_f[threadIdx.x] = static_cast<float>(rho_weight[j_load]);
        omega_tile_f[threadIdx.x] = static_cast<float>(omega[j_load]);
        kappa_tile_f[threadIdx.x] = static_cast<float>(kappa[j_load]);
      }
    }
    __syncthreads();

    if(active_rho_i) {
      const auto tile_npts = min(blockDim.x, npts - tile);
      if constexpr(UseFloatPair) {
        for(int j = 0; j < tile_npts; ++j) {
          const auto r_ij2 = distance2(coords_tile[j], r_i);
          const auto r_ij2_f = static_cast<float>(r_ij2);
          const auto g_ij = fmaf(r_ij2_f, omega_i_f, kappa_i_f);
          const auto g_ji = fmaf(r_ij2_f, omega_tile_f[j], kappa_tile_f[j]);
          const auto g_sum = g_ij + g_ji;
          const auto g_ij_inv = __fdividef(1.0f, g_ij);
          const auto g_ji_inv = __fdividef(1.0f, g_ji);
          const auto g_sum_inv = __fdividef(1.0f, g_sum);
          auto pair = rho_weight_tile_f[j] * g_ij_inv * g_ji_inv * g_sum_inv;
          F_i += static_cast<double>(pair);
          pair *= g_ij_inv + g_sum_inv;
          U_i += static_cast<double>(pair);
          W_i += static_cast<double>(pair * r_ij2_f);
        }
      } else {
        for(int j = 0; j < tile_npts; ++j) {
          const auto r_ij2 = distance2(coords_tile[j], r_i);
          const auto g_ij = r_ij2 * omega_i + kappa_i;
          const auto g_ji = r_ij2 * omega_tile[j] + kappa_tile[j];
          const auto g_sum = g_ij + g_ji;
          const auto g_ij_inv = 1.0 / g_ij;
          const auto g_ji_inv = 1.0 / g_ji;
          const auto g_sum_inv = 1.0 / g_sum;
          auto pair = rho_weight_tile[j] * g_ij_inv * g_ji_inv * g_sum_inv;
          F_i += pair;
          pair *= g_ij_inv + g_sum_inv;
          U_i += pair;
          W_i += pair * r_ij2;
        }
      }
    }
    __syncthreads();
  }

  if(active_i) {
    if(!active_rho_i) {
      eps[i] = 0.0;
      vrho[i] = 0.0;
      vgamma[i] = 0.0;
      return;
    }

    const auto gamma_i = gamma[i];
    const auto rho_1 = 1.0 / rho_i;
    const auto rho_2 = rho_1 * rho_1;
    const auto rho_4 = rho_2 * rho_2;
    const auto rho_5 = rho_1 * rho_4;
    const auto gamma2 = gamma_i * gamma_i;
    const auto omega2 = params.c * gamma2 * rho_4 + four_pi_over_three * rho_i;
    const auto omega_eval = sqrt(omega2);
    const auto domega_drho = 0.5 * (four_pi_over_three - 4.0 * params.c * gamma2 * rho_5) / omega_eval;
    const auto domega_dgamma = params.c * gamma_i * rho_4 / omega_eval;
    const auto scaled_dkappa_drho = kappa_i / 6.0;
    const auto F = -1.5 * F_i;

    eps[i] = cuda_beta(params) + 0.5 * F;
    vrho[i] = cuda_beta(params) + F + 1.5 * (U_i * scaled_dkappa_drho + W_i * rho_i * domega_drho);
    vgamma[i] = 1.5 * W_i * rho_i * domega_dgamma;
  }
}

__global__ void grid_gradient_kernel(int npts, const double* coords,
  const double* omega, const double* kappa, int vv_npts,
  const double* vv_coords, const double* vv_rho_weight,
  const double* vv_omega, const double* vv_kappa, double* gradient) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= npts) return;

  const auto r_i = load_point(coords, i);
  double gx = 0.0;
  double gy = 0.0;
  double gz = 0.0;

  for(int j = 0; j < vv_npts; ++j) {
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

__global__ void grid_gradient_excluding_same_parent_kernel(int npts,
  const double* coords, const double* omega, const double* kappa,
  const double* rho_weight, const std::int32_t* parent, double* gradient) {

  __shared__ double3 coords_tile[threads_per_block];
  __shared__ double rho_weight_tile[threads_per_block];
  __shared__ double omega_tile[threads_per_block];
  __shared__ double kappa_tile[threads_per_block];
  __shared__ std::int32_t parent_tile[threads_per_block];

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  const auto active_i = i < npts;
  const auto r_i = active_i ? load_point(coords, i) : double3{0.0, 0.0, 0.0};
  const auto omega_i = active_i ? omega[i] : 0.0;
  const auto kappa_i = active_i ? kappa[i] : 0.0;
  const auto parent_i = active_i ? parent[i] : std::int32_t{-1};

  double gx = 0.0;
  double gy = 0.0;
  double gz = 0.0;

  for(int tile = 0; tile < npts; tile += blockDim.x) {
    const auto j_load = tile + threadIdx.x;
    if(j_load < npts) {
      coords_tile[threadIdx.x] = load_point(coords, j_load);
      rho_weight_tile[threadIdx.x] = rho_weight[j_load];
      omega_tile[threadIdx.x] = omega[j_load];
      kappa_tile[threadIdx.x] = kappa[j_load];
      parent_tile[threadIdx.x] = parent[j_load];
    }
    __syncthreads();

    if(active_i) {
      const auto tile_npts = min(blockDim.x, npts - tile);
      for(int j = 0; j < tile_npts; ++j) {
        if(parent_i == parent_tile[j]) continue;
        const auto r_j = coords_tile[j];
        const auto dx = r_j.x - r_i.x;
        const auto dy = r_j.y - r_i.y;
        const auto dz = r_j.z - r_i.z;
        const auto r_ij2 = dx*dx + dy*dy + dz*dz;
        const auto g_ij = r_ij2 * omega_i + kappa_i;
        const auto g_ji = r_ij2 * omega_tile[j] + kappa_tile[j];
        const auto g_sum = g_ij + g_ji;
        const auto g_ij_inv = 1.0 / g_ij;
        const auto g_ji_inv = 1.0 / g_ji;
        const auto g_sum_inv = 1.0 / g_sum;
        const auto pair = rho_weight_tile[j] * g_ij_inv * g_ji_inv * g_sum_inv;
        const auto prefactor = pair * (omega_i * g_ij_inv + omega_tile[j] * g_ji_inv
          + (omega_i + omega_tile[j]) * g_sum_inv);
        gx += prefactor * dx;
        gy += prefactor * dy;
        gz += prefactor * dz;
      }
    }
    __syncthreads();
  }

  if(active_i) {
    gradient[3*i + 0] = -3.0 * gx;
    gradient[3*i + 1] = -3.0 * gy;
    gradient[3*i + 2] = -3.0 * gz;
  }
}

__global__ void hessian_intermediates_kernel(int npts, const double* coords,
  const double* weights, const double* rho, const double* omega,
  const double* kappa, double* U, double* W, double* A, double* B,
  double* C, double* E) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= npts) return;

  const auto r_i = load_point(coords, i);
  double U_i = 0.0;
  double W_i = 0.0;
  double A_i = 0.0;
  double B_i = 0.0;
  double C_i = 0.0;
  double E_i = 0.0;
  for(int j = 0; j < npts; ++j) {
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

__global__ void hessian_contraction_kernel(int npts, int ntrial,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa, const double* U,
  const double* W, const double* A, const double* B, const double* C,
  const double* domega_drho, const double* domega_dgamma,
  const double* dkappa_drho, const double* d2omega_drho2,
  const double* d2omega_dgamma2, const double* d2omega_drho_dgamma,
  const double* d2kappa_drho2, const double* rho_t,
  const double* gamma_t, double* f_rho_t, double* f_gamma_t) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  const auto trial = blockIdx.y;
  if(i >= npts || trial >= ntrial) return;

  const auto r_i = load_point(coords, i);
  double f_rho_i = 0.0;
  double f_gamma_i = 0.0;
  for(int j = 0; j < npts; ++j) {
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
    const auto offset_j = trial * npts + j;
    f_rho_i += weights[j] * (f_rho_rho_ij * rho_t[offset_j] + f_rho_gamma_ij * gamma_t[offset_j]);
    f_gamma_i += weights[j] * (f_gamma_rho_ij * rho_t[offset_j] + f_gamma_gamma_ij * gamma_t[offset_j]);
  }

  const auto f_rho_rho_ii = 2.0 * domega_drho[i] * W[i] + 2.0 * dkappa_drho[i] * U[i]
    + rho[i] * (d2omega_drho2[i] * W[i] + d2kappa_drho2[i] * U[i]
    + dkappa_drho[i] * dkappa_drho[i] * A[i] + domega_drho[i] * domega_drho[i] * C[i]
    + 2.0 * domega_drho[i] * dkappa_drho[i] * B[i]);
  const auto f_gamma_rho_ii = domega_dgamma[i] * W[i] + rho[i] * (d2omega_drho_dgamma[i] * W[i]
    + domega_dgamma[i] * (dkappa_drho[i] * B[i] + domega_drho[i] * C[i]));
  const auto f_gamma_gamma_ii = rho[i] * (d2omega_dgamma2[i] * W[i]
    + domega_dgamma[i] * domega_dgamma[i] * C[i]);
  const auto offset_i = trial * npts + i;
  f_rho_i += f_rho_rho_ii * rho_t[offset_i] + f_gamma_rho_ii * gamma_t[offset_i];
  f_gamma_i += f_gamma_rho_ii * rho_t[offset_i] + f_gamma_gamma_ii * gamma_t[offset_i];
  f_rho_t[offset_i] = f_rho_i;
  f_gamma_t[offset_i] = f_gamma_i;
}

__global__ void grid_response_kernel(int npts, int natoms,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa,
  const int* grid_associated_atom, double* Egr, double* Ugr, double* Wgr) {

  const auto i = blockIdx.x * blockDim.x + threadIdx.x;
  const auto atom = blockIdx.y;
  if(i >= npts || atom >= natoms) return;
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
    return;
  }

  const auto r_i = load_point(coords, i);
  const auto i_in_atom = i_atom == atom;
  double3 Eg{0.0, 0.0, 0.0};
  double3 Ug{0.0, 0.0, 0.0};
  double3 Wg{0.0, 0.0, 0.0};
  for(int j = 0; j < npts; ++j) {
    const auto j_atom = grid_associated_atom[j];
    if(j_atom < 0) continue;
    const auto j_in_atom = j_atom == atom;
    if((!i_in_atom && !j_in_atom) || (i_in_atom && j_in_atom)) continue;
    const auto r_j = load_point(coords, j);
    const double3 r_ji{r_j.x - r_i.x, r_j.y - r_i.y, r_j.z - r_i.z};
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

void check_last_error(const char* msg) {
  const auto err = cudaGetLastError();
  GAUXC_CUDA_ERROR(msg, err);
}

}

void eval_omega_kappa_cuda(cudaStream_t stream, std::int32_t npts,
  const double* rho, const double* gamma, Parameters params, double* omega,
  double* kappa, double* domega_drho, double* domega_dgamma,
  double* d2omega_drho2, double* d2omega_dgamma2,
  double* d2omega_drho_dgamma, double* dkappa_drho,
  double* d2kappa_drho2) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  omega_kappa_kernel<<<blocks, threads_per_block, 0, stream>>>(npts, rho, gamma,
    params, omega, kappa, domega_drho, domega_dgamma, d2omega_drho2,
    d2omega_dgamma2, d2omega_drho_dgamma, dkappa_drho, d2kappa_drho2);
  check_last_error("VV10 omega/kappa kernel failed");
}

void eval_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* omega, const double* kappa,
  std::int32_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* F, double* U, double* W) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  fock_kernel<<<blocks, threads_per_block, 0, stream>>>(npts, coords, omega,
    kappa, vv_npts, vv_coords, vv_rho_weight, vv_omega, vv_kappa, F, U, W);
  check_last_error("VV10 fock kernel failed");
}

void eval_exc_vxc_from_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* rho, const double* gamma, const double* F, const double* U,
  const double* W, Parameters params, double* eps, double* vrho,
  double* vgamma) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  exc_vxc_from_fock_kernel<<<blocks, threads_per_block, 0, stream>>>(npts, rho,
    gamma, F, U, W, params, eps, vrho, vgamma);
  check_last_error("VV10 EXC/VXC finalize kernel failed");
}

void eval_exc_vxc_fock_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* rho, const double* gamma,
  const double* omega, const double* kappa, const double* rho_weight,
  Parameters params, double* eps, double* vrho, double* vgamma) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  if(params.pair_math_mode == 1) {
    exc_vxc_fock_kernel<true><<<blocks, threads_per_block, 0, stream>>>(npts, coords,
      rho, gamma, omega, kappa, rho_weight, params, eps, vrho, vgamma);
  } else {
    exc_vxc_fock_kernel<false><<<blocks, threads_per_block, 0, stream>>>(npts, coords,
      rho, gamma, omega, kappa, rho_weight, params, eps, vrho, vgamma);
  }
  check_last_error("VV10 fused EXC/VXC fock kernel failed");
}

void eval_grid_gradient_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* omega, const double* kappa,
  std::int32_t vv_npts, const double* vv_coords,
  const double* vv_rho_weight, const double* vv_omega,
  const double* vv_kappa, double* gradient) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  grid_gradient_kernel<<<blocks, threads_per_block, 0, stream>>>(npts, coords,
    omega, kappa, vv_npts, vv_coords, vv_rho_weight, vv_omega, vv_kappa,
    gradient);
  check_last_error("VV10 grid gradient kernel failed");
}

void eval_grid_gradient_excluding_same_parent_cuda(cudaStream_t stream,
  std::int32_t npts, const double* coords, const double* omega,
  const double* kappa, const double* rho_weight,
  const std::int32_t* parent, double* gradient) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  grid_gradient_excluding_same_parent_kernel<<<blocks, threads_per_block, 0, stream>>>(
    npts, coords, omega, kappa, rho_weight, parent, gradient);
  check_last_error("VV10 parent-excluding grid gradient kernel failed");
}

void eval_hessian_intermediates_cuda(cudaStream_t stream, std::int32_t npts,
  const double* coords, const double* weights, const double* rho,
  const double* omega, const double* kappa, double* U, double* W,
  double* A, double* B, double* C, double* E) {
  const auto blocks = util::div_ceil(npts, threads_per_block);
  hessian_intermediates_kernel<<<blocks, threads_per_block, 0, stream>>>(npts,
    coords, weights, rho, omega, kappa, U, W, A, B, C, E);
  check_last_error("VV10 hessian intermediates kernel failed");
}

void eval_hessian_contraction_cuda(cudaStream_t stream, std::int32_t npts,
  std::int32_t ntrial, const double* coords, const double* weights,
  const double* rho, const double* omega, const double* kappa,
  const double* U, const double* W, const double* A, const double* B,
  const double* C, const double* domega_drho,
  const double* domega_dgamma, const double* dkappa_drho,
  const double* d2omega_drho2, const double* d2omega_dgamma2,
  const double* d2omega_drho_dgamma, const double* d2kappa_drho2,
  const double* rho_t, const double* gamma_t, double* f_rho_t,
  double* f_gamma_t) {
  const dim3 blocks(util::div_ceil(npts, threads_per_block), ntrial);
  hessian_contraction_kernel<<<blocks, threads_per_block, 0, stream>>>(npts,
    ntrial, coords, weights, rho, omega, kappa, U, W, A, B, C,
    domega_drho, domega_dgamma, dkappa_drho, d2omega_drho2,
    d2omega_dgamma2, d2omega_drho_dgamma, d2kappa_drho2, rho_t,
    gamma_t, f_rho_t, f_gamma_t);
  check_last_error("VV10 hessian contraction kernel failed");
}

void eval_grid_response_cuda(cudaStream_t stream, std::int32_t npts,
  std::int32_t natoms, const double* coords, const double* weights,
  const double* rho, const double* omega, const double* kappa,
  const std::int32_t* grid_associated_atom, double* Egr, double* Ugr,
  double* Wgr) {
  const dim3 blocks(util::div_ceil(npts, threads_per_block), natoms);
  grid_response_kernel<<<blocks, threads_per_block, 0, stream>>>(npts, natoms,
    coords, weights, rho, omega, kappa, grid_associated_atom, Egr, Ugr, Wgr);
  check_last_error("VV10 grid response kernel failed");
}

}

#endif