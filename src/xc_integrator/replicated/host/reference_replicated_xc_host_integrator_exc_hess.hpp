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

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"
#include "host/gauxc_hess_kernel.hpp"
#include "host/util.hpp"
#include "common/integrator_constants.hpp"
#include <gauxc/molmeta.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace GauXC::detail {

/**
 *  Log-derivatives of one grid point's partition weight w = q P_C / Z, in
 *  REDUCED coordinates: every atom but the parent C, the point held fixed.
 *  The parent's rows and columns follow by translational invariance.
 *
 *    dlw  = (d w) / w  = a_C - g,        g = sum_D pi_D a_D
 *    d2lw = (d2 w) / w = (a_C-g)(a_C-g)^T + B_C
 *                        - sum_D pi_D [ (a_D-g)(a_D-g)^T + B_D ]
 *
 *  with pi_D = P_D / Z, a_D = d ln P_D = sum_E t_DE dmu_DE and
 *  B_D = sum_E (u_DE - t_DE^2) dmu_DE dmu_DE^T + t_DE d2mu_DE, where
 *  t = s'/s and u = s''/s of the cell function. Every term stays finite
 *  as P_C -> 0. mu_DE with its derivatives, and s, t, u, are GENERATED
 *  (xckernel gauxcwriter) -- t and u from the factored ln s, which does
 *  not cancel as mu -> 1.
 *
 *  Returns false when the derivatives vanish identically: an SSF point
 *  inside the parent's cutoff sphere, or a zero partition.
 */
struct WeightDerivScratch {
  std::vector<double> rd, P, S, T, U, a, g;
};

inline bool partition_weight_log_derivs( bool is_becke, const Molecule& mol,
  const MolMeta& meta, int32_t iC, double dist_nearest, const double* rg,
  WeightDerivScratch& w, double* dlw, double* d2lw ) {

  const int32_t N  = static_cast<int32_t>(mol.natoms());
  const int32_t n3 = 3*N;
  const auto&  RAB = meta.rab();
  const double a_ssf = integrator::magic_ssf_factor<>;
  w.rd.resize(N); w.P.resize(N);
  w.S.resize(size_t(N)*N); w.T.resize(size_t(N)*N); w.U.resize(size_t(N)*N);
  w.a.resize(size_t(N)*n3); w.g.resize(n3);

  for( int32_t D = 0; D < N; ++D ) {
    const double dx = rg[0]-mol[D].x, dy = rg[1]-mol[D].y, dz = rg[2]-mol[D].z;
    w.rd[D] = std::sqrt( dx*dx + dy*dy + dz*dz );
  }
  // SSF: inside this sphere the partition is exactly 1, flat in every atom
  if( not is_becke and w.rd[iC] < 0.5*(1.-a_ssf)*dist_nearest ) return false;

  for( int32_t D = 0; D < N; ++D ) {
    w.P[D] = 1.;
    for( int32_t E = 0; E < N; ++E ) {
      if( E == D ) continue;
      const double mu = (w.rd[D] - w.rd[E]) / RAB[E + D*N];
      double s = 1., t = 0., u = 0.;
      if( is_becke ? (1. - mu < 1e-10) : (mu >= a_ssf) ) {
        s = 0.;                        // flat zero (Becke: the point sits on E)
      } else if( is_becke or mu > -a_ssf ) {
        #include "host/gauxc_hess_call_cell.inc"
      }
      w.S[size_t(D)*N+E] = s; w.T[size_t(D)*N+E] = t; w.U[size_t(D)*N+E] = u;
      w.P[D] *= s;
    }
  }
  double Z = 0.;
  for( int32_t D = 0; D < N; ++D ) Z += w.P[D];
  if( Z <= 0. or w.P[iC] <= 0. ) return false;

  std::fill( w.g.begin(), w.g.end(), 0. );
  std::fill( d2lw, d2lw + size_t(n3)*n3, 0. );
  auto active = [&]( int32_t D ) { return D == iC or w.P[D] > 1e-16*Z; };

  for( int32_t D = 0; D < N; ++D ) {
    if( not active(D) ) continue;
    double* aD = w.a.data() + size_t(D)*n3;
    std::fill( aD, aD + n3, 0. );
    const double piD   = w.P[D] / Z;
    const double bcoef = (D == iC ? 1. : 0.) - piD;
    for( int32_t E = 0; E < N; ++E ) {
      if( E == D ) continue;
      const double t = w.T[size_t(D)*N+E], u = w.U[size_t(D)*N+E];
      if( t == 0. and u == 0. ) continue;
      const double RD[3] = { mol[D].x, mol[D].y, mol[D].z };
      const double RE[3] = { mol[E].x, mol[E].y, mol[E].z };
      double mu, dmu[6], d2mu[36];
      #include "host/gauxc_hess_call_mu.inc"
      (void)mu;
      int idx[6];
      for( int c = 0; c < 3; ++c ) {
        idx[c]   = (D == iC) ? -1 : 3*D + c;   // the parent is not a free
        idx[3+c] = (E == iC) ? -1 : 3*E + c;   // coordinate here
      }
      for( int i = 0; i < 6; ++i ) if( idx[i] >= 0 ) aD[idx[i]] += t * dmu[i];
      if( bcoef != 0. )
      for( int i = 0; i < 6; ++i ) if( idx[i] >= 0 )
      for( int j = 0; j < 6; ++j ) if( idx[j] >= 0 )
        d2lw[ size_t(idx[i])*n3 + idx[j] ] +=
          bcoef * ( (u - t*t)*dmu[i]*dmu[j] + t*d2mu[6*i+j] );
    }
    for( int32_t k = 0; k < n3; ++k ) w.g[k] += piD * aD[k];
  }

  const double* aC = w.a.data() + size_t(iC)*n3;
  for( int32_t k = 0; k < n3; ++k ) dlw[k] = aC[k] - w.g[k];
  for( int32_t k = 0; k < n3; ++k )
  for( int32_t l = 0; l < n3; ++l ) d2lw[ size_t(k)*n3 + l ] += dlw[k]*dlw[l];
  for( int32_t D = 0; D < N; ++D ) {
    if( not active(D) ) continue;
    const double  piD = w.P[D] / Z;
    const double* aD  = w.a.data() + size_t(D)*n3;
    for( int32_t k = 0; k < n3; ++k ) {
      const double vk = piD * (aD[k] - w.g[k]);
      if( vk == 0. ) continue;
      for( int32_t l = 0; l < n3; ++l )
        d2lw[ size_t(k)*n3 + l ] -= vk * (aD[l] - w.g[l]);
    }
  }
  return true;
}

/**
 *  RKS nuclear Hessian of the XC energy at fixed density matrix.
 *
 *  The XC energy depends on the nuclear coordinates three ways: the basis
 *  functions ride their atoms, the grid points ride their parent atoms, and
 *  the partition weights change. IntegratorSettingsEXC_HESS mirrors
 *  exc_grad:
 *
 *   - include_weight_derivatives = false: the basis class alone -- fixed
 *     grid, fixed weights, the analogue of exc_grad's Hellmann-Feynman
 *     option.
 *   - include_weight_derivatives = true (default): all three. A task's
 *     points ride its parent C, so its energy depends on R_A - R_C only;
 *     each task is accumulated in coordinates that exclude C and C's rows
 *     and columns are restored by translational invariance. The weight
 *     class adds w'' e + w' e'^T + e' w'^T per point, with the Becke or
 *     SSF weight derivatives from partition_weight_log_derivs (LKO is not
 *     implemented).
 *
 *  The per-point kernels are GENERATED (xckernel/emitters/gauxcwriter.py from
 *  engine/geometric.geometric_hessian), so the contraction cannot drift from
 *  the expressions that tool validates against finite differences. The
 *  assembly recipe here -- sum the per-function rows over the shells of one
 *  atom to get the nuclear-perturbed fields, contract the two row sets as an
 *  OUTER PRODUCT, then add the Pulay and same-atom seeds -- is itself checked
 *  against the unfactorised expression by that tool, because every kernel can
 *  be individually correct while the recipe is a wrong reading of them.
 */
template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_hess_( int64_t m, int64_t n, const value_type* P, int64_t ldp,
                  value_type* EXC_HESS, const IntegratorSettingsXC& ks_settings ) {

  const auto& basis = this->load_balancer_->basis();
  const int64_t nbf = basis.nbf();
  if( m != n )   GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf ) GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");
  if( ldp < nbf ) GAUXC_GENERIC_EXCEPTION("Invalid LDP");

  this->load_balancer_->get_tasks();

  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_hess_local_work_( P, ldp, EXC_HESS, ks_settings );
  });

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");
    const int n3 = 3*static_cast<int>(this->load_balancer_->molecule().natoms());
    this->reduction_driver_->allreduce_inplace( EXC_HESS, n3*n3, ReductionOp::Sum );
  });
}


template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_hess_local_work_( const value_type* P, int64_t ldp,
                        value_type* EXC_HESS, const IntegratorSettingsXC& settings ) {

  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();
  const auto& func  = *this->func_;

  if( func.needs_laplacian() )
    GAUXC_GENERIC_EXCEPTION("EXC Hessian Not Implemented For Laplacian-Dependent MGGAs");
  if( func.is_polarized() )
    GAUXC_GENERIC_EXCEPTION("EXC Hessian Only Implemented For RKS");

  BasisSetMap basis_map(basis, mol);
  const int32_t nbf    = basis.nbf();
  const int32_t natoms = static_cast<int32_t>(mol.natoms());
  const int32_t n3     = 3*natoms;

  auto& tasks = this->load_balancer_->get_tasks();

  for( int32_t i = 0; i < n3*n3; ++i ) EXC_HESS[i] = 0.;

  const bool is_gga  = func.is_gga() or func.is_mgga();
  const bool is_mgga = func.is_mgga();
  const int  nfield  = is_mgga ? 3 : (is_gga ? 2 : 1); // rho[, sigma[, tau]]

  IntegratorSettingsEXC_HESS hess_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsEXC_HESS*>(&settings) )
    hess_settings = *tmp;
  const bool full = hess_settings.include_weight_derivatives;

  const auto& molmeta = this->load_balancer_->molmeta();
  const auto& lb_state = this->load_balancer_->state();
  if( full and not lb_state.modified_weights_are_stored )
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");
  const XCWeightAlg weight_alg = lb_state.weight_alg;
  if( full and weight_alg == XCWeightAlg::LKO )
    GAUXC_GENERIC_EXCEPTION("EXC Hessian Weight Derivatives Not Implemented For LKO");
  const bool is_becke    = weight_alg == XCWeightAlg::Becke;
  const bool partitioned = weight_alg != XCWeightAlg::NOTPARTITIONED;

  const size_t ntasks = tasks.size();
  #pragma omp parallel
  {
  XCHostData<value_type> host_data;
  std::vector<value_type> hess_local(size_t(n3)*n3, 0.);
  // this task's contribution, before the translational-invariance fold
  std::vector<value_type> HR(size_t(n3)*n3);
  std::vector<double> dlw(n3), d2lw( full ? size_t(n3)*n3 : 0 );
  WeightDerivScratch wscr;

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    auto& task = tasks[iT];
    std::fill( HR.begin(), HR.end(), 0. );

    const int32_t npts    = static_cast<int32_t>(task.points.size());
    const int32_t nbe     = task.bfn_screening.nbe;
    const int32_t nshells = static_cast<int32_t>(task.bfn_screening.shell_list.size());
    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // ---- scratch ------------------------------------------------------
    // basis + grad(3) + hess(6) + der3(10): the Hessian needs one more
    // derivative than the gradient does, because d2 chi / dA dB of the
    // GRADIENT ingredient is a third derivative of the basis function.
    // Even an LDA needs the basis HESSIAN: the delta_AB term carries a
    // double displacement of one function. A GGA needs der3 on top,
    // because d2/dAdB of the GRADIENT ingredient is a third derivative.
    const int ncomp = is_gga ? 20 : 10;
    host_data.basis_eval.resize( size_t(ncomp) * npts * nbe );
    host_data.nbe_scr   .resize( size_t(nbe) * nbe );
    host_data.zmat      .resize( size_t(4) * npts * nbe );
    host_data.den_scr   .resize( size_t(4) * npts );
    host_data.eps       .resize( npts );
    host_data.vrho      .resize( npts );
    host_data.v2rho2    .resize( npts );
    if( is_gga ) {
      host_data.gamma     .resize( npts );
      host_data.vgamma    .resize( npts );
      host_data.v2rhogamma.resize( npts );
      host_data.v2gamma2  .resize( npts );
    }
    if( is_mgga ) {
      host_data.tau       .resize( npts );
      host_data.vtau      .resize( npts );
      host_data.v2rhotau  .resize( npts );
      host_data.v2gammatau.resize( npts );
      host_data.v2tau2    .resize( npts );
    }

    auto* basis_eval = host_data.basis_eval.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* xmat       = host_data.zmat.data();
    auto* den_eval   = host_data.den_scr.data();

    auto* dbf_x = basis_eval + size_t(1)*npts*nbe;
    auto* dbf_y = basis_eval + size_t(2)*npts*nbe;
    auto* dbf_z = basis_eval + size_t(3)*npts*nbe;
    value_type* d2[6] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
    value_type* d3[10];
    for( int i = 0; i < 10; ++i ) d3[i] = nullptr;
    for( int i = 0; i < 6; ++i ) d2[i] = basis_eval + size_t(4+i)*npts*nbe;
    if( is_gga )
      for( int i = 0; i < 10; ++i ) d3[i] = basis_eval + size_t(10+i)*npts*nbe;

    auto [submat_map, foo] =
      gen_compressed_submat_map( basis_map, task.bfn_screening.shell_list, nbf, nbf );

    if( is_gga )
      lwd->eval_collocation_der3( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbf_x, dbf_y, dbf_z, d2[0],d2[1],d2[2],d2[3],d2[4],d2[5],
        d3[0],d3[1],d3[2],d3[3],d3[4],d3[5],d3[6],d3[7],d3[8],d3[9] );
    else
      lwd->eval_collocation_hessian( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbf_x, dbf_y, dbf_z, d2[0],d2[1],d2[2],d2[3],d2[4],d2[5] );

    // X = 2 P B and its gradient rows; the 4 blocks are contiguous, matching
    // the contiguous basis+gradient layout eval_xmat expects.
    const double xmat_fac = 2.0;
    lwd->eval_xmat( 4*npts, nbf, nbe, submat_map, xmat_fac, P, ldp, basis_eval, nbe,
                    xmat, nbe, nbe_scr );
    auto* xmat_x = xmat + size_t(1)*npts*nbe;
    auto* xmat_y = xmat + size_t(2)*npts*nbe;
    auto* xmat_z = xmat + size_t(3)*npts*nbe;

    auto* dden_x = den_eval + size_t(1)*npts;
    auto* dden_y = den_eval + size_t(2)*npts;
    auto* dden_z = den_eval + size_t(3)*npts;
    auto* gamma  = host_data.gamma.data();
    auto* tau    = host_data.tau.data();

    if( is_mgga )
      lwd->eval_uvvar_mgga_rks( npts, nbe, basis_eval, dbf_x, dbf_y, dbf_z, nullptr,
        xmat, nbe, xmat_x, xmat_y, xmat_z, nbe, den_eval, dden_x, dden_y, dden_z,
        gamma, tau, nullptr );
    else if( is_gga )
      lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbf_x, dbf_y, dbf_z,
        xmat, nbe, den_eval, dden_x, dden_y, dden_z, gamma );
    else
      lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, xmat, nbe, den_eval );

    auto* vrho   = host_data.vrho.data();
    auto* vgamma = host_data.vgamma.data();
    auto* vtau   = host_data.vtau.data();
    auto* v2rho2      = host_data.v2rho2.data();
    auto* v2rhogamma  = host_data.v2rhogamma.data();
    auto* v2gamma2    = host_data.v2gamma2.data();
    auto* v2rhotau    = host_data.v2rhotau.data();
    auto* v2gammatau  = host_data.v2gammatau.data();
    auto* v2tau2      = host_data.v2tau2.data();

    if( is_mgga )
      func.eval_vxc_fxc( npts, den_eval, gamma, nullptr, tau, vrho, vgamma,
        nullptr, vtau, v2rho2, v2rhogamma, nullptr, v2rhotau, v2gamma2,
        nullptr, v2gammatau, nullptr, nullptr, v2tau2 );
    else if( is_gga )
      func.eval_vxc_fxc( npts, den_eval, gamma, vrho, vgamma,
                         v2rho2, v2rhogamma, v2gamma2 );
    else
      func.eval_vxc_fxc( npts, den_eval, vrho, v2rho2 );

    // the energy density itself, for the weight-class terms
    auto* eps = host_data.eps.data();
    if( full and partitioned ) {
      if( is_mgga )     func.eval_exc( npts, den_eval, gamma, nullptr, tau, eps );
      else if( is_gga ) func.eval_exc( npts, den_eval, gamma, eps );
      else              func.eval_exc( npts, den_eval, eps );
    }

    // ---- atoms touched by this task -----------------------------------
    std::vector<int32_t> atoms;               // local index -> global atom
    std::vector<int32_t> sh_atom(nshells);    // shell -> LOCAL atom index
    for( int32_t ish = 0; ish < nshells; ++ish ) {
      const int iAt = basis_map.shell_to_center( shell_list[ish] );
      auto it = std::find( atoms.begin(), atoms.end(), iAt );
      if( it == atoms.end() ) { sh_atom[ish] = static_cast<int32_t>(atoms.size());
                                atoms.push_back(iAt); }
      else                    { sh_atom[ish] = static_cast<int32_t>(it - atoms.begin()); }
    }
    const int32_t nat_loc = static_cast<int32_t>(atoms.size());

    // rows[(atom,dir)][field][pt]: the nuclear-perturbed fields. Summing the
    // generated per-function rows over the shells of one atom is what turns
    // the non-Pulay pair term into a rank update.
    const int nrow = nfield + 3;              // F_* plus G_x,G_y,G_z
    std::vector<double> rows( size_t(nat_loc)*3*nrow*npts, 0. );
    auto ROW = [&]( int a, int d, int k, int ip ) -> double& {
      return rows[ (((size_t(a)*3 + d)*nrow + k)*npts) + ip ];
    };

    size_t bf_off = 0;
    for( int32_t ish = 0; ish < nshells; ++ish ) {
      const int sh_sz = basis[shell_list[ish]].size();
      const int a     = sh_atom[ish];
      for( int ibf = 0, mu = static_cast<int>(bf_off); ibf < sh_sz; ++ibf, ++mu ) {
        for( int32_t ip = 0; ip < npts; ++ip ) {
          const size_t k = size_t(mu) + size_t(ip)*nbe;
          const double U0 = xmat[k];
          const double U1 = is_gga ? xmat_x[k] : 0.0;
          const double U2 = is_gga ? xmat_y[k] : 0.0;
          const double U3 = is_gga ? xmat_z[k] : 0.0;
          const double db[3] = { dbf_x[k], dbf_y[k], dbf_z[k] };
          // d2[] packing is xx,xy,xz,yy,yz,zz
          const int h_idx[3][3] = { {0,1,2}, {1,3,4}, {2,4,5} };
          for( int d = 0; d < 3; ++d ) {
            // the -d/dr sign of a displaced basis function, folded in
            const double dchi = -db[d];
            double ddchi[3] = {0.,0.,0.};
            if( is_gga )
              for( int i = 0; i < 3; ++i ) ddchi[i] = -d2[ h_idx[d][i] ][k];
            double F_rho=0., F_sigma=0., F_tau=0., Gx=0., Gy=0., Gz=0.;
            // generated: kernel call AND the row writes, so the row
            // layout here and in the pair reads below come from one table
            #include "host/gauxc_hess_call_rows.inc"
          }
        }
      }
      bf_off += sh_sz;
    }

    // ---- outer-product (non-Pulay) pair term --------------------------
    for( int32_t a = 0; a < nat_loc; ++a )
    for( int32_t b = 0; b < nat_loc; ++b )
    for( int dx = 0; dx < 3; ++dx )
    for( int dy = 0; dy < 3; ++dy ) {
      double acc = 0.;
      for( int32_t ip = 0; ip < npts; ++ip ) {
        double h = 0.;
        #include "host/gauxc_hess_call_pair.inc"
        acc += weights[ip] * h;
      }
      HR[ size_t(3*atoms[a]+dx)*n3 + (3*atoms[b]+dy) ] += acc;
    }

    const int h_idx[3][3] = { {0,1,2}, {1,3,4}, {2,4,5} };

    // ---- same-atom (delta_AB) term -------------------------------------
    // Both displacements hit the SAME function: no pair, no factorisation,
    // but only nbe * 9 * npts work.
    {
      // d3 packing from eval_collocation_der3: xxx,xxy,xxz,xyy,xyz,xzz,
      //                                        yyy,yyz,yzz,zzz
      static const int d3_idx[3][3][3] = {{{0,1,2},{1,3,4},{2,4,5}},
                                          {{1,3,4},{3,6,7},{4,7,8}},
                                          {{2,4,5},{4,7,8},{5,8,9}}};
      size_t off_u = 0;
      for( int32_t ish = 0; ish < nshells; ++ish ) {
        const int sh_u = basis[shell_list[ish]].size();
        const int a    = sh_atom[ish];
        for( int iu = 0, mu = static_cast<int>(off_u); iu < sh_u; ++iu, ++mu )
        for( int dx = 0; dx < 3; ++dx )
        for( int dy = 0; dy < 3; ++dy ) {
          double acc = 0.;
          for( int32_t ip = 0; ip < npts; ++ip ) {
            const size_t k = size_t(mu) + size_t(ip)*nbe;
            // two displacement signs cancel: this carries +d2/dxdy
            const double d2c = d2[ h_idx[dx][dy] ][k];
            double d3c[3] = {0.,0.,0.};
            if( is_gga )
              for( int i = 0; i < 3; ++i ) d3c[i] = d3[ d3_idx[dx][dy][i] ][k];
            double sv = 0.;
            #include "host/gauxc_hess_call_same.inc"
            acc += sv;
          }
          HR[ size_t(3*atoms[a]+dx)*n3 + (3*atoms[a]+dy) ] += acc;
        }
        off_u += sh_u;
      }
    }

    // ---- Pulay term, as matrix products --------------------------------
    // The generated weights W_st(g) make the Pulay term bilinear in the
    // displaced rows of the two functions,
    //   pulay_uv = sum_g sum_st Phi^dx_s(u,g) W_st(g) Phi^dy_t(v,g),
    // slot 0 the displaced function, slots 1-3 its displaced gradient. So
    // M^{dx,dy} = Phi^dx (W Phi^dy)^T is one GEMM, and the Hessian block
    // is its D_uv-weighted sum over the functions of each atom pair.
    {
      std::vector<double> Ploc( size_t(nbe)*nbe );
      detail::submat_set( nbf, nbf, nbe, nbe, P, ldp, Ploc.data(), nbe, submat_map );

      std::vector<int32_t> bf_atom( nbe );
      for( int32_t ish = 0, mu = 0; ish < nshells; ++ish )
        for( int i = 0; i < basis[shell_list[ish]].size(); ++i )
          bf_atom[mu++] = atoms[ sh_atom[ish] ];

      // per-point weights, Wbuf[(4*s+t)*npts + ip]
      std::vector<double> Wbuf( size_t(16)*npts, 0. );
      for( int32_t ip = 0; ip < npts; ++ip ) {
        double W[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        #include "host/gauxc_hess_call_pulayW.inc"
        for( int i = 0; i < 16; ++i ) Wbuf[ size_t(i)*npts + ip ] = W[i];
      }

      // displaced rows, the -d/dr sign folded in: nbe x (ns*npts) per direction
      const int ns = is_gga ? 4 : 1;
      const size_t blk = size_t(npts)*nbe;
      std::vector<double> Phi( 3*ns*blk ), Z( ns*blk ), M( size_t(nbe)*nbe );
      auto PHI = [&]( int d, int s ) { return Phi.data() + (size_t(d)*ns + s)*blk; };
      value_type* dbf[3] = { dbf_x, dbf_y, dbf_z };
      for( int d = 0; d < 3; ++d )
      for( size_t k = 0; k < blk; ++k ) {
        PHI(d,0)[k] = -dbf[d][k];
        for( int i = 0; i < ns-1; ++i ) PHI(d,i+1)[k] = -d2[ h_idx[d][i] ][k];
      }

      for( int dy = 0; dy < 3; ++dy ) {
        // Z_s = sum_t W_st Phi^dy_t, column (point) scaling
        std::fill( Z.begin(), Z.end(), 0. );
        for( int s = 0; s < ns; ++s )
        for( int t = 0; t < ns; ++t ) {
          const double* w  = Wbuf.data() + size_t(4*s+t)*npts;
          const double* ph = PHI(dy,t);
          double*       z  = Z.data() + s*blk;
          for( int32_t ip = 0; ip < npts; ++ip ) {
            if( w[ip] == 0. ) continue;
            for( int32_t mu = 0; mu < nbe; ++mu )
              z[ mu + size_t(ip)*nbe ] += w[ip] * ph[ mu + size_t(ip)*nbe ];
          }
        }
        for( int dx = 0; dx < 3; ++dx ) {
          blas::gemm( 'N', 'T', nbe, nbe, ns*npts, 1., PHI(dx,0), nbe,
                      Z.data(), nbe, 0., M.data(), nbe );
          for( int32_t nu = 0; nu < nbe; ++nu )
          for( int32_t mu = 0; mu < nbe; ++mu )
            HR[ size_t(3*bf_atom[mu]+dx)*n3 + (3*bf_atom[nu]+dy) ] +=
              xmat_fac * Ploc[ size_t(mu) + size_t(nu)*nbe ] * M[ size_t(mu) + size_t(nu)*nbe ];
        }
      }
    }

    // ---- weight class ------------------------------------------------
    //   d2(w e) = w'' e + w' e'^T + e' w'^T + w e''
    // in reduced coordinates; w e'' is the basis class above. e' is the
    // basis-class energy-density gradient, from the same per-atom rows.
    if( full and partitioned ) {
      const int32_t iC = task.iParent;
      std::vector<double> de_loc( size_t(nat_loc)*3*npts );
      for( int32_t a = 0; a < nat_loc; ++a )
      for( int d = 0; d < 3; ++d )
      for( int32_t ip = 0; ip < npts; ++ip ) {
        double de = 0.;
        #include "host/gauxc_hess_call_egrad.inc"
        de_loc[ (size_t(a)*3 + d)*npts + ip ] = de;
      }
      for( int32_t ip = 0; ip < npts; ++ip ) {
        if( weights[ip] == 0. ) continue;
        if( not partition_weight_log_derivs( is_becke, mol, molmeta, iC,
              task.dist_nearest, points + 3*size_t(ip), wscr, dlw.data(), d2lw.data() ) )
          continue;
        const double we = weights[ip] * den_eval[ip] * eps[ip];
        for( size_t k = 0; k < size_t(n3)*n3; ++k ) HR[k] += we * d2lw[k];
        for( int32_t a = 0; a < nat_loc; ++a ) {
          if( atoms[a] == iC ) continue;       // C's functions ride the grid
          for( int d = 0; d < 3; ++d ) {
            const double v = weights[ip] * de_loc[ (size_t(a)*3 + d)*npts + ip ];
            const size_t K = 3*size_t(atoms[a]) + d;
            for( int32_t J = 0; J < n3; ++J ) {
              HR[ size_t(J)*n3 + K ] += dlw[J] * v;
              HR[ K*n3 + J ]         += v * dlw[J];
            }
          }
        }
      }
    }

    // ---- fold into the Hessian -----------------------------------------
    // Basis-only: as is. Full: the task's points ride the parent C, so its
    // energy depends on R_A - R_C only. Drop C's rows and columns (C's own
    // functions ride the grid) and restore them by translational invariance.
    if( not full ) {
      for( size_t k = 0; k < size_t(n3)*n3; ++k ) hess_local[k] += HR[k];
    } else {
      const int32_t iC = task.iParent;
      for( int32_t A = 0; A < natoms; ++A ) if( A != iC )
      for( int dx = 0; dx < 3; ++dx )
      for( int32_t B = 0; B < natoms; ++B ) if( B != iC )
      for( int dy = 0; dy < 3; ++dy ) {
        const double v = HR[ size_t(3*A+dx)*n3 + (3*B+dy) ];
        if( v == 0. ) continue;
        hess_local[ size_t(3*A +dx)*n3 + (3*B +dy) ] += v;
        hess_local[ size_t(3*A +dx)*n3 + (3*iC+dy) ] -= v;
        hess_local[ size_t(3*iC+dx)*n3 + (3*B +dy) ] -= v;
        hess_local[ size_t(3*iC+dx)*n3 + (3*iC+dy) ] += v;
      }
    }

  } // tasks

  #pragma omp critical
  {
    for( int32_t i = 0; i < n3*n3; ++i ) EXC_HESS[i] += hess_local[i];
  }
  } // omp parallel
}

} // namespace GauXC::detail
