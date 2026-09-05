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
#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/sycl_util.hpp"
#include "device/xc_device_data.hpp"
#include "sycl_extensions.hpp"
#include "sycl_atomics.hpp"
#include "sycl_launch.hpp"

#define MGGA_KERNEL_SM_BLOCK 32

namespace GauXC {



// den_shared has 3+!!need_lapl rows, each [warp_size][MGGA_KERNEL_SM_BLOCK+1].
// The row count depends on the need_lapl template parameter; since need_lapl
// is a template (compile-time) parameter here, the extent "3+!!need_lapl" is
// itself a compile-time constant and local_mem() can be sized from it
// directly, exactly like the CUDA
// __shared__ double den_shared[3+!!need_lapl][warp_size][MGGA_KERNEL_SM_BLOCK+1]
// declaration.
template <bool trial, density_id den_select, bool need_lapl>
void eval_vvar_mgga_kern( size_t           ntasks,
                           XCDeviceTask* tasks_device) {

  auto it = GauXC::sycl::this_item();
  constexpr auto warp_size = sycl::warp_size;
  auto& den_shared = GauXC::sycl::local_mem<double[3+(need_lapl?1:0)][warp_size][MGGA_KERNEL_SM_BLOCK+1]>();

  const int batch_idx = it.get_group(0);
  if( batch_idx >= ntasks ) return;

  auto& task = tasks_device[ batch_idx ];

  const auto npts            = task.npts;
  const auto nbf             = task.bfn_screening.nbe;
  double* tau_eval_device  = nullptr;
  double* lapl_eval_device = nullptr;

  if constexpr (trial){
    if constexpr (den_select == DEN_S) {
      tau_eval_device = task.ttau_s;
      if constexpr (need_lapl) {
        lapl_eval_device = task.tlapl_s;
      }
    }
    if constexpr (den_select == DEN_Z) {
      tau_eval_device = task.ttau_z;
      if constexpr (need_lapl) {
        lapl_eval_device = task.tlapl_z;
      }
    }
  } else{
    if constexpr (den_select == DEN_S) {
      tau_eval_device = task.tau_s;
      if constexpr (need_lapl) {
        lapl_eval_device = task.lapl_s;
      }
    }
    if constexpr (den_select == DEN_Z) {
      tau_eval_device = task.tau_z;
      if constexpr (need_lapl) {
        lapl_eval_device = task.lapl_z;
      }
    }
  }

  //const auto* basis_eval_device = task.bf;
  const auto* dbasis_x_eval_device = task.dbfx;
  const auto* dbasis_y_eval_device = task.dbfy;
  const auto* dbasis_z_eval_device = task.dbfz;
  decltype(dbasis_x_eval_device) basis_lapl_eval_device = nullptr;
  if constexpr (need_lapl) {
    basis_lapl_eval_device = task.d2bflapl;
  }

  //const auto* den_basis_prod_device    = task.zmat;
  const auto* den_basis_dx_prod_device = task.xmat_x;
  const auto* den_basis_dy_prod_device = task.xmat_y;
  const auto* den_basis_dz_prod_device = task.xmat_z;
  decltype(den_basis_dx_prod_device) den_basis_prod_device = nullptr;
  if constexpr (need_lapl) {
    den_basis_prod_device = task.zmat;
  }

  const int threadIdx_x = it.get_local_id(2);
  const int threadIdx_y = it.get_local_id(1);
  const int blockDim_x  = it.get_local_range(2);
  const int blockDim_y  = it.get_local_range(1);

  for ( int bid_x = it.get_group(2) * blockDim_x;
        bid_x < nbf;
        bid_x += blockDim_x * it.get_group_range(2) ) {

    for ( int bid_y = it.get_group(1) * MGGA_KERNEL_SM_BLOCK;
          bid_y < npts;
          bid_y += MGGA_KERNEL_SM_BLOCK * it.get_group_range(1) ) {

      for (int sm_y = threadIdx_y; sm_y < MGGA_KERNEL_SM_BLOCK; sm_y += blockDim_y) {
        den_shared[0][threadIdx_x][sm_y] = 0.;
        den_shared[1][threadIdx_x][sm_y] = 0.;
        den_shared[2][threadIdx_x][sm_y] = 0.;
        if constexpr (need_lapl)
          den_shared[3][threadIdx_x][sm_y] = 0.;

        if (bid_y + threadIdx_x < npts and bid_x + sm_y < nbf) {
          const double* db_x_col = den_basis_dx_prod_device + (bid_x + sm_y)*npts;
          const double* db_y_col = den_basis_dy_prod_device + (bid_x + sm_y)*npts;
          const double* db_z_col = den_basis_dz_prod_device + (bid_x + sm_y)*npts;

          const double* bf_x_col = dbasis_x_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_y_col = dbasis_y_eval_device  + (bid_x + sm_y)*npts;
          const double* bf_z_col = dbasis_z_eval_device  + (bid_x + sm_y)*npts;


          den_shared[0][threadIdx_x][sm_y] = bf_x_col[ bid_y + threadIdx_x ] * db_x_col[ bid_y + threadIdx_x ];
          den_shared[1][threadIdx_x][sm_y] = bf_y_col[ bid_y + threadIdx_x ] * db_y_col[ bid_y + threadIdx_x ];
          den_shared[2][threadIdx_x][sm_y] = bf_z_col[ bid_y + threadIdx_x ] * db_z_col[ bid_y + threadIdx_x ];


          if constexpr (need_lapl) {
            const double* db_col   = den_basis_prod_device  + (bid_x + sm_y)*npts;
            const double* bf_l_col = basis_lapl_eval_device + (bid_x + sm_y)*npts;
            den_shared[3][threadIdx_x][sm_y] = bf_l_col[ bid_y + threadIdx_x ] * db_col[ bid_y + threadIdx_x ];
          }
        }
      }
      ::sycl::group_barrier( it.get_group() );


      for (int sm_y = threadIdx_y; sm_y < MGGA_KERNEL_SM_BLOCK; sm_y += blockDim_y) {
        const int tid_y = bid_y + sm_y;

        double tx_reg  = den_shared[0][sm_y][threadIdx_x];
        double ty_reg  = den_shared[1][sm_y][threadIdx_x];
        double tz_reg  = den_shared[2][sm_y][threadIdx_x];
        // Warp blocks are stored col major
        double tau_reg = 0.0;
        tau_reg  = 0.5 * sycl::warp_reduce_sum<warp_size>( tx_reg );
        tau_reg += 0.5 * sycl::warp_reduce_sum<warp_size>( ty_reg );
        tau_reg += 0.5 * sycl::warp_reduce_sum<warp_size>( tz_reg );

        double lapl_reg = 0.0;
        if constexpr (need_lapl) {
          lapl_reg = den_shared[3][sm_y][threadIdx_x];
          lapl_reg = sycl::warp_reduce_sum<warp_size>( lapl_reg );
          lapl_reg = 2. * lapl_reg + 4. * tau_reg;
        }

        if( threadIdx_x == 0 and tid_y < npts ) {
          atomic_add_device( tau_eval_device   + tid_y, tau_reg );
          if constexpr (need_lapl) {
            atomic_add_device( lapl_eval_device   + tid_y, lapl_reg );
          }
        }
      }
      ::sycl::group_barrier( it.get_group() );
    }
  }
}




template <bool need_lapl>
void eval_uvars_mgga_uks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {

  auto it = GauXC::sycl::this_item();
  const int batch_idx = it.get_group(0);
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts  = task.npts;

  auto* tau_pos_eval_device = task.tau_s;
  auto* tau_neg_eval_device = task.tau_z;

  double* lapl_pos_eval_device = nullptr;
  double* lapl_neg_eval_device = nullptr;
  if constexpr (need_lapl) {
    lapl_pos_eval_device = task.lapl_s;
    lapl_neg_eval_device = task.lapl_z;
  }

  const int tid = it.get_global_id(2);

  if( tid < npts ) {
    const double ts = tau_pos_eval_device[ tid ];
    const double tz = tau_neg_eval_device[ tid ];
    tau_pos_eval_device[ tid ] = 0.5*(ts + tz);
    tau_neg_eval_device[ tid ] = 0.5*(ts - tz);

    if constexpr (need_lapl) {
      const double ls = lapl_pos_eval_device[ tid ];
      const double lz = lapl_neg_eval_device[ tid ];
      lapl_pos_eval_device[ tid ] = 0.5*(ls + lz);
      lapl_neg_eval_device[ tid ] = 0.5*(ls - lz);
    }
  }

}


inline void eval_tmat_mgga_rks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {
  auto it = GauXC::sycl::this_item();
  const int batch_idx = it.get_group(0);
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts  = task.npts;

  const auto*   dden_sx_eval_device = task.dden_sx;
  const auto*   dden_sy_eval_device = task.dden_sy;
  const auto*   dden_sz_eval_device = task.dden_sz;
  const auto*   tdden_sx_eval_device = task.tdden_sx;
  const auto*   tdden_sy_eval_device = task.tdden_sy;
  const auto*   tdden_sz_eval_device = task.tdden_sz;

  const auto* weight_device  = task.weights;
  const auto* vgamma_device  = task.vgamma;
  const auto* v2rho2_device     = task.v2rho2;
  const auto* v2rhogamma_device = task.v2rhogamma;
  const auto* v2gamma2_device   = task.v2gamma2;
  const auto* v2rhotau_device  = task.v2rhotau;
  const auto* v2tau2_device   = task.v2tau2;
  const auto* v2gammatau_device = task.v2gammatau;
  const auto* trho_device       = task.tden_s;
  const auto* ttau_device       = task.ttau_s;

  auto* FXC_A_device   = task.FXC_A_s;
  auto* FXC_Bx_device   = task.FXC_Bx_s;
  auto* FXC_By_device   = task.FXC_By_s;
  auto* FXC_Bz_device   = task.FXC_Bz_s;
  auto* FXC_C_device   = task.FXC_C_s;

  const int tid = it.get_global_id(2);

  if( tid < npts ) {
    const auto dx = dden_sx_eval_device[ tid ];
    const auto dy = dden_sy_eval_device[ tid ];
    const auto dz = dden_sz_eval_device[ tid ];
    const auto tdx = tdden_sx_eval_device[ tid ];
    const auto tdy = tdden_sy_eval_device[ tid ];
    const auto tdz = tdden_sz_eval_device[ tid ];
    const auto tgamma = tdx*dx + tdy*dy + tdz*dz;

    const auto FXC_A = v2rho2_device[ tid ] * trho_device[ tid ] + 2.0 * v2rhogamma_device[tid] * tgamma +
        v2rhotau_device[ tid ] * ttau_device[ tid ];
    FXC_A_device[ tid ]  = weight_device[ tid ] * FXC_A;

    const auto FXC_C = v2rhotau_device[ tid ] * trho_device[ tid ] + 2.0 * v2gammatau_device[ tid ] * tgamma +
        v2tau2_device[ tid ] * ttau_device[ tid ];
    FXC_C_device[ tid ]  = weight_device[ tid ] * FXC_C;

    const auto B_coef = v2rhogamma_device[tid] * trho_device[tid] + 2.0 * v2gamma2_device[tid] * tgamma +
        v2gammatau_device[ tid ] * ttau_device[ tid ];
    FXC_Bx_device[ tid ] = 2.0 * weight_device[ tid ] * ( B_coef * dx + vgamma_device[ tid ] * tdx );
    FXC_By_device[ tid ] = 2.0 * weight_device[ tid ] * ( B_coef * dy + vgamma_device[ tid ] * tdy );
    FXC_Bz_device[ tid ] = 2.0 * weight_device[ tid ] * ( B_coef * dz + vgamma_device[ tid ] * tdz );
  }

}



inline void eval_tmat_mgga_uks_kernel( size_t ntasks, XCDeviceTask* tasks_device) {

  auto it = GauXC::sycl::this_item();
  const int batch_idx = it.get_group(0);
  if( batch_idx >= ntasks ) return;

  const auto& task = tasks_device[ batch_idx ];
  const auto npts            = task.npts;

  const auto* tden_s_device   = task.tden_s;
  const auto* tden_z_device   = task.tden_z;
  const auto* ttau_s_device   = task.ttau_s;
  const auto* ttau_z_device   = task.ttau_z;
  const auto* weight_device   = task.weights;

  const auto*     tden_pos_x_eval_device = task.tdden_sx;
  const auto*     tden_pos_y_eval_device = task.tdden_sy;
  const auto*     tden_pos_z_eval_device = task.tdden_sz;
  const auto*     den_pos_x_eval_device = task.dden_sx;
  const auto*     den_pos_y_eval_device = task.dden_sy;
  const auto*     den_pos_z_eval_device = task.dden_sz;

  const auto*     tden_neg_x_eval_device = task.tdden_zx;
  const auto*     tden_neg_y_eval_device = task.tdden_zy;
  const auto*     tden_neg_z_eval_device = task.tdden_zz;
  const auto*     den_neg_x_eval_device = task.dden_zx;
  const auto*     den_neg_y_eval_device = task.dden_zy;
  const auto*     den_neg_z_eval_device = task.dden_zz;

  const double* vgamma_aa_device   = task.vgamma_pp;
  const double* vgamma_ab_device   = task.vgamma_pm;
  const double* vgamma_bb_device   = task.vgamma_mm;
  const double* v2rho2_a_a_device    = task.v2rho2_a_a;
  const double* v2rho2_a_b_device    = task.v2rho2_a_b;
  const double* v2rho2_b_b_device    = task.v2rho2_b_b;
  const double* v2rhogamma_a_aa_device = task.v2rhogamma_a_aa;
  const double* v2rhogamma_a_ab_device = task.v2rhogamma_a_ab;
  const double* v2rhogamma_a_bb_device = task.v2rhogamma_a_bb;
  const double* v2rhogamma_b_aa_device = task.v2rhogamma_b_aa;
  const double* v2rhogamma_b_ab_device = task.v2rhogamma_b_ab;
  const double* v2rhogamma_b_bb_device = task.v2rhogamma_b_bb;
  const double* v2gamma2_aa_aa_device = task.v2gamma2_aa_aa;
  const double* v2gamma2_aa_ab_device = task.v2gamma2_aa_ab;
  const double* v2gamma2_aa_bb_device = task.v2gamma2_aa_bb;
  const double* v2gamma2_ab_ab_device = task.v2gamma2_ab_ab;
  const double* v2gamma2_ab_bb_device = task.v2gamma2_ab_bb;
  const double* v2gamma2_bb_bb_device = task.v2gamma2_bb_bb;
  const double* v2rhotau_a_a_device   = task.v2rhotau_a_a;
  const double* v2rhotau_a_b_device   = task.v2rhotau_a_b;
  const double* v2rhotau_b_a_device   = task.v2rhotau_b_a;
  const double* v2rhotau_b_b_device   = task.v2rhotau_b_b;
  const double* v2gammatau_aa_a_device= task.v2gammatau_aa_a;
  const double* v2gammatau_aa_b_device= task.v2gammatau_aa_b;
  const double* v2gammatau_ab_a_device= task.v2gammatau_ab_a;
  const double* v2gammatau_ab_b_device= task.v2gammatau_ab_b;
  const double* v2gammatau_bb_a_device= task.v2gammatau_bb_a;
  const double* v2gammatau_bb_b_device= task.v2gammatau_bb_b;
  const double* v2tau2_a_a_device   = task.v2tau2_a_a;
  const double* v2tau2_a_b_device   = task.v2tau2_a_b;
  const double* v2tau2_b_b_device   = task.v2tau2_b_b;

  auto* FXC_A_s_device        = task.FXC_A_s;
  auto* FXC_A_z_device        = task.FXC_A_z;
  auto* FXC_Bx_s_device       = task.FXC_Bx_s;
  auto* FXC_Bx_z_device       = task.FXC_Bx_z;
  auto* FXC_By_s_device       = task.FXC_By_s;
  auto* FXC_By_z_device       = task.FXC_By_z;
  auto* FXC_Bz_s_device       = task.FXC_Bz_s;
  auto* FXC_Bz_z_device       = task.FXC_Bz_z;
  auto* FXC_C_s_device        = task.FXC_C_s;
  auto* FXC_C_z_device        = task.FXC_C_z;

  const int tid = it.get_global_id(2);

  if( tid < npts ) {
    const auto ps = tden_s_device[ tid ];
    const auto pz = tden_z_device[ tid ];
    const auto trho_a_device = 0.5*(ps + pz);
    const auto trho_b_device = 0.5*(ps - pz);
    const auto ts = ttau_s_device[ tid ];
    const auto tz = ttau_z_device[ tid ];
    const auto tau_a = 0.5*(ts + tz);
    const auto tau_b = 0.5*(ts - tz);

    const auto tdndx   = tden_pos_x_eval_device[ tid ];
    const auto tdndy   = tden_pos_y_eval_device[ tid ];
    const auto tdndz   = tden_pos_z_eval_device[ tid ];
    const auto tdMzdx  = tden_neg_x_eval_device[ tid ];
    const auto tdMzdy  = tden_neg_y_eval_device[ tid ];
    const auto tdMzdz  = tden_neg_z_eval_device[ tid ];
    const auto tdden_a_x = 0.5*(tdndx + tdMzdx);
    const auto tdden_a_y = 0.5*(tdndy + tdMzdy);
    const auto tdden_a_z = 0.5*(tdndz + tdMzdz);
    const auto tdden_b_x = 0.5*(tdndx - tdMzdx);
    const auto tdden_b_y = 0.5*(tdndy - tdMzdy);
    const auto tdden_b_z = 0.5*(tdndz - tdMzdz);

    const auto dndx   = den_pos_x_eval_device[ tid ];
    const auto dndy   = den_pos_y_eval_device[ tid ];
    const auto dndz   = den_pos_z_eval_device[ tid ];
    const auto dMzdx  = den_neg_x_eval_device[ tid ];
    const auto dMzdy  = den_neg_y_eval_device[ tid ];
    const auto dMzdz  = den_neg_z_eval_device[ tid ];
    const auto dden_a_x = 0.5*(dndx + dMzdx);
    const auto dden_a_y = 0.5*(dndy + dMzdy);
    const auto dden_a_z = 0.5*(dndz + dMzdz);
    const auto dden_b_x = 0.5*(dndx - dMzdx);
    const auto dden_b_y = 0.5*(dndy - dMzdy);
    const auto dden_b_z = 0.5*(dndz - dMzdz);

    const auto tgamma_pp = tdden_a_x * dden_a_x + tdden_a_y * dden_a_y + tdden_a_z * dden_a_z;
    const auto tgamma_pm = tdden_a_x * dden_b_x + tdden_a_y * dden_b_y + tdden_a_z * dden_b_z
                                 + tdden_b_x * dden_a_x + tdden_b_y * dden_a_y + tdden_b_z * dden_a_z;
    const auto tgamma_mm = tdden_b_x * dden_b_x + tdden_b_y * dden_b_y + tdden_b_z * dden_b_z;


    const auto A_a = v2rho2_a_a_device[tid] * trho_a_device + 2.0 * v2rhogamma_a_aa_device[tid] * tgamma_pp +
          v2rhogamma_a_ab_device[tid] * tgamma_pm + 2.0 * v2rhogamma_a_bb_device[tid] * tgamma_mm +
          v2rho2_a_b_device[tid] * trho_b_device + v2rhotau_a_a_device[tid] * tau_a +
          v2rhotau_a_b_device[tid] * tau_b;
    const auto A_b = v2rho2_b_b_device[tid] * trho_b_device + 2.0 * v2rhogamma_b_bb_device[tid] * tgamma_mm +
          v2rhogamma_b_ab_device[tid] * tgamma_pm + 2.0 * v2rhogamma_b_aa_device[tid] * tgamma_pp +
          v2rho2_a_b_device[tid] * trho_a_device + v2rhotau_b_b_device[tid] * tau_b +
          v2rhotau_b_a_device[tid] * tau_a;
    FXC_A_s_device[ tid ] = 0.5 * weight_device[ tid ] * (A_a + A_b);
    FXC_A_z_device[ tid ] = 0.5 * weight_device[ tid ] * (A_a - A_b);

    // Compute C coefficients for alpha and beta spin
    const auto C_a = v2rhotau_a_a_device[tid] * trho_a_device + v2rhotau_b_a_device[tid] * trho_b_device
             + 2.0 * v2gammatau_aa_a_device[tid] * tgamma_pp  + v2gammatau_ab_a_device[tid] * tgamma_pm
             + 2.0 * v2gammatau_bb_a_device[tid] * tgamma_mm
             + v2tau2_a_a_device[tid] * tau_a + v2tau2_a_b_device[tid] * tau_b;

    const auto C_b = v2rhotau_a_b_device[tid] * trho_a_device + v2rhotau_b_b_device[tid] * trho_b_device
             + 2.0 * v2gammatau_aa_b_device[tid] * tgamma_pp + v2gammatau_ab_b_device[tid] * tgamma_pm
             + 2.0 * v2gammatau_bb_b_device[tid] * tgamma_mm
             + v2tau2_a_b_device[tid] * tau_a + v2tau2_b_b_device[tid] * tau_b;

    FXC_C_s_device[tid] = 0.5 * weight_device[tid] * (C_a + C_b);
    FXC_C_z_device[tid] = 0.5 * weight_device[tid] * (C_a - C_b);

    // Calculate B coefficients for alpha spin
    const double B_coef1_a = v2rhogamma_a_aa_device[tid] * trho_a_device   + 2.0 * v2gamma2_aa_aa_device[tid] * tgamma_pp +
                 v2gamma2_aa_ab_device[tid] * tgamma_pm + 2.0 * v2gamma2_aa_bb_device[tid] * tgamma_mm +
                 v2rhogamma_b_aa_device[tid] * trho_b_device + v2gammatau_aa_a_device[tid] * tau_a +
                 v2gammatau_aa_b_device[tid] * tau_b;

    const double B_coef2_a = v2rhogamma_a_ab_device[tid] * trho_a_device + 2.0 * v2gamma2_aa_ab_device[tid] * tgamma_pp +
          v2gamma2_ab_ab_device[tid] * tgamma_pm + 2.0 * v2gamma2_ab_bb_device[tid] * tgamma_mm +
          v2rhogamma_b_ab_device[tid] * trho_b_device + v2gammatau_ab_a_device[tid] * tau_a +
          v2gammatau_ab_b_device[tid] * tau_b;

    // Calculate gradient components for alpha spin
    const double Bx_a = 2.0 * B_coef1_a * dden_a_x + B_coef2_a * dden_b_x +
           2.0 * vgamma_aa_device[tid] * tdden_a_x + vgamma_ab_device[tid] * tdden_b_x;

    const double By_a = 2.0 * B_coef1_a * dden_a_y + B_coef2_a * dden_b_y +
           2.0 * vgamma_aa_device[tid] * tdden_a_y + vgamma_ab_device[tid] * tdden_b_y;

    const double Bz_a = 2.0 * B_coef1_a * dden_a_z + B_coef2_a * dden_b_z +
           2.0 * vgamma_aa_device[tid] * tdden_a_z + vgamma_ab_device[tid] * tdden_b_z;

    // Calculate B coefficients for beta spin
    const double B_coef1_b = v2rhogamma_b_bb_device[tid] * trho_b_device + 2.0 * v2gamma2_bb_bb_device[tid] * tgamma_mm +
          v2gamma2_ab_bb_device[tid] * tgamma_pm + 2.0 * v2gamma2_aa_bb_device[tid] * tgamma_pp +
          v2rhogamma_a_bb_device[tid] * trho_a_device + v2gammatau_bb_b_device[tid] * tau_b +
          v2gammatau_bb_a_device[tid] * tau_a;

    const double B_coef2_b = v2rhogamma_b_ab_device[tid] * trho_b_device + 2.0 * v2gamma2_ab_bb_device[tid] * tgamma_mm +
          v2gamma2_ab_ab_device[tid] * tgamma_pm + 2.0 * v2gamma2_aa_ab_device[tid] * tgamma_pp +
          v2rhogamma_a_ab_device[tid] * trho_a_device + v2gammatau_ab_b_device[tid] * tau_b +
          v2gammatau_ab_a_device[tid] * tau_a;

    const double Bx_b = 2.0 * B_coef1_b * dden_b_x + B_coef2_b * dden_a_x +
           2.0 * vgamma_bb_device[tid] * tdden_b_x + vgamma_ab_device[tid] * tdden_a_x;

    const double By_b = 2.0 * B_coef1_b * dden_b_y + B_coef2_b * dden_a_y +
           2.0 * vgamma_bb_device[tid] * tdden_b_y + vgamma_ab_device[tid] * tdden_a_y;

    const double Bz_b = 2.0 * B_coef1_b * dden_b_z + B_coef2_b * dden_a_z +
           2.0 * vgamma_bb_device[tid] * tdden_b_z + vgamma_ab_device[tid] * tdden_a_z;

    // Store weighted values in output arrays
    FXC_Bx_s_device[tid] = 0.5 * weight_device[tid] * (Bx_a + Bx_b);
    FXC_By_s_device[tid] = 0.5 * weight_device[tid] * (By_a + By_b);
    FXC_Bz_s_device[tid] = 0.5 * weight_device[tid] * (Bz_a + Bz_b);
    FXC_Bx_z_device[tid] = 0.5 * weight_device[tid] * (Bx_a - Bx_b);
    FXC_By_z_device[tid] = 0.5 * weight_device[tid] * (By_a - By_b);
    FXC_Bz_z_device[tid] = 0.5 * weight_device[tid] * (Bz_a - Bz_b);

  }

}


}
