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
#include "host/gauxc_nc_kernel.hpp"
#include <cmath>
#include <vector>

namespace GauXC::detail {

/**
 *  GKS (two-component, noncollinear) FXC contraction.
 *
 *  The energy is the locally collinear one exc_vxc evaluates for GKS:
 *  n_+- = (rho_s +- |m|)/2 fed to the spin-polarized functional. The
 *  kernel applied to a trial density, per noncollinear field slot, is
 *  GENERATED (xckernel ncwriter: the mechanical second derivative of
 *  that map), and assembled exactly as the potential is -- so
 *  FXC_X = d/dh VXC_X(P + h tP) for X = s, z, y, x.
 *
 *  Below gks_dtol the map is not twice differentiable (the transverse
 *  kernel carries 1/|m|); there the generated collinear limit is used:
 *  every magnetization component responds like the spin channel of a
 *  collinear perturbation about the spin-symmetric reference.
 *
 *  LDA only for now. The GGA map of Scalmani and Frisch has a second
 *  singularity, 1/|(grad rho_s . grad m_J)_J|, which vanishes at density
 *  critical points at finite |m|; its regularization is still open.
 *
 *  Argument order follows the GKS eval_exc_vxc: (s, z, y, x).
 */
template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_fxc_contraction_( int64_t m, int64_t n,
                         const value_type* Ps, int64_t ldps, const value_type* Pz, int64_t ldpz,
                         const value_type* Py, int64_t ldpy, const value_type* Px, int64_t ldpx,
                         const value_type* tPs, int64_t ldtps, const value_type* tPz, int64_t ldtpz,
                         const value_type* tPy, int64_t ldtpy, const value_type* tPx, int64_t ldtpx,
                         value_type* FXCs, int64_t ldfxcs, value_type* FXCz, int64_t ldfxcz,
                         value_type* FXCy, int64_t ldfxcy, value_type* FXCx, int64_t ldfxcx,
                         const IntegratorSettingsXC& ks_settings ) {

  const auto& basis = this->load_balancer_->basis();
  const int64_t nbf = basis.nbf();
  if( m != n )   GAUXC_GENERIC_EXCEPTION("P/FXC Must Be Square");
  if( m != nbf ) GAUXC_GENERIC_EXCEPTION("P/FXC Must Have Same Dimension as Basis");
  for( int64_t ld : { ldps, ldpz, ldpy, ldpx, ldtps, ldtpz, ldtpy, ldtpx,
                      ldfxcs, ldfxcz, ldfxcy, ldfxcx } )
    if( ld < nbf ) GAUXC_GENERIC_EXCEPTION("Invalid Leading Dimension");

  // Symmetrize the trial densities: Exc sees only the symmetric part of a
  // density matrix, but the gradient channel of a noncollinear GGA would
  // not (see #225); done here so the contraction never depends on it.
  std::vector<value_type> tsym[4];
  const value_type* tPin[4] = { tPs, tPz, tPy, tPx };
  const int64_t     ldtin[4] = { ldtps, ldtpz, ldtpy, ldtpx };
  const value_type* tP[4];
  int64_t ldtP[4];
  for( int k = 0; k < 4; ++k ) {
    tsym[k].resize( nbf*nbf );
    for( int64_t j = 0; j < nbf; ++j )
      for( int64_t i = 0; i < nbf; ++i )
        tsym[k][i + j*nbf] = 0.5 * ( tPin[k][i + j*ldtin[k]] + tPin[k][j + i*ldtin[k]] );
    tP[k] = tsym[k].data(); ldtP[k] = nbf;
  }
  const value_type* P[4]   = { Ps, Pz, Py, Px };
  const int64_t     ldP[4] = { ldps, ldpz, ldpy, ldpx };
  value_type* FXC[4]       = { FXCs, FXCz, FXCy, FXCx };
  const int64_t ldFXC[4]   = { ldfxcs, ldfxcz, ldfxcy, ldfxcx };

  auto& tasks = this->load_balancer_->get_tasks();
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    fxc_contraction_gks_local_work_( basis, P, ldP, tP, ldtP, FXC, ldFXC,
                                     ks_settings, tasks.begin(), tasks.end() );
  });

  this->timer_.time_op("XCIntegrator.Allreduce", [&](){
    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");
    for( int k = 0; k < 4; ++k )
      this->reduction_driver_->allreduce_inplace( FXC[k], nbf*nbf, ReductionOp::Sum );
  });
}


template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  fxc_contraction_gks_local_work_( const basis_type& basis,
                                   const value_type* const P[4], const int64_t ldP[4],
                                   const value_type* const tP[4], const int64_t ldtP[4],
                                   value_type* const FXC[4], const int64_t ldFXC[4],
                                   const IntegratorSettingsXC& settings,
                                   task_iterator task_begin, task_iterator task_end ) {

  IntegratorSettingsKS ks_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsKS*>(&settings) )
    ks_settings = *tmp;
  const double dtol = ks_settings.gks_dtol;

  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());
  const auto& func = *this->func_;
  const auto& mol  = this->load_balancer_->molecule();

  if( not func.is_polarized() )
    GAUXC_GENERIC_EXCEPTION("GKS FXC Contraction Requires A Polarized Functional");
  if( not func.is_lda() )
    GAUXC_GENERIC_EXCEPTION("GKS FXC Contraction Only Implemented For LDA");

  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored )
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Been Modified");

  BasisSetMap basis_map(basis, mol);
  const int32_t nbf = basis.nbf();

  for( int k = 0; k < 4; ++k )
    for( int32_t j = 0; j < nbf; ++j )
      for( int32_t i = 0; i < nbf; ++i )
        FXC[k][i + j*ldFXC[k]] = 0.;

  const size_t ntasks = std::distance(task_begin, task_end);

  #pragma omp parallel
  {
  std::vector<value_type> basis_eval, nbe_scr, X, tX, Z, fields, tfields,
                          npm, vrho, v2rho2;

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {
    const auto& task = *(task_begin + iT);
    const int32_t npts    = static_cast<int32_t>(task.points.size());
    const int32_t nbe     = task.bfn_screening.nbe;
    const int32_t nshells = static_cast<int32_t>(task.bfn_screening.shell_list.size());
    const auto* points    = task.points.data()->data();
    const auto* weights   = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();
    const size_t blk = size_t(npts)*nbe;

    basis_eval.resize( blk );  nbe_scr.resize( size_t(nbe)*nbe );
    X.resize( 4*blk ); tX.resize( 4*blk ); Z.resize( 4*blk );
    fields.resize( 4*size_t(npts) ); tfields.resize( 4*size_t(npts) );
    npm.resize( 2*size_t(npts) ); vrho.resize( 2*size_t(npts) ); v2rho2.resize( 3*size_t(npts) );

    std::vector< std::array<int32_t,3> > submat_map;
    std::tie(submat_map, std::ignore) =
      gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list, basis_eval.data() );

    // fields and trial fields, slots (s, z, y, x); GKS uses xmat_fac = 1
    for( int k = 0; k < 4; ++k ) {
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, P[k], ldP[k], basis_eval.data(), nbe,
                      X.data() + k*blk, nbe, nbe_scr.data() );
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, tP[k], ldtP[k], basis_eval.data(), nbe,
                      tX.data() + k*blk, nbe, nbe_scr.data() );
      for( int32_t ip = 0; ip < npts; ++ip ) {
        const auto* b = basis_eval.data() + size_t(ip)*nbe;
        fields [k*npts + ip] = blas::dot( nbe, b, 1, X .data() + k*blk + size_t(ip)*nbe, 1 );
        tfields[k*npts + ip] = blas::dot( nbe, b, 1, tX.data() + k*blk + size_t(ip)*nbe, 1 );
      }
    }

    // the locally collinear variables at the TRUE |m| (not the vxc fallback's
    // component average): the limit kernel below needs n+ = n- at m = 0
    for( int32_t ip = 0; ip < npts; ++ip ) {
      const double mz = fields[1*npts+ip], my = fields[2*npts+ip], mx = fields[3*npts+ip];
      const double mn = std::sqrt( mx*mx + my*my + mz*mz );
      npm[2*ip]   = 0.5*( fields[ip] + mn );
      npm[2*ip+1] = 0.5*( fields[ip] - mn );
    }
    func.eval_vxc_fxc( npts, npm.data(), vrho.data(), v2rho2.data() );

    // kernel per point, in the potential's field slots; Z_X = 1/2 k_X chi
    for( int32_t ip = 0; ip < npts; ++ip ) {
      const double rs = fields[ip], rz = fields[npts+ip], ry = fields[2*npts+ip], rx = fields[3*npts+ip];
      const double ts = tfields[ip], tz = tfields[npts+ip], ty = tfields[2*npts+ip], tx = tfields[3*npts+ip];
      const double v0 = vrho[2*ip], v1 = vrho[2*ip+1];
      const double f0 = v2rho2[3*ip], f1 = v2rho2[3*ip+1], f2 = v2rho2[3*ip+2];
      double ks, kx, ky, kz;
      if( std::sqrt( rx*rx + ry*ry + rz*rz ) > dtol )
        xckernel::nc_fxc_contract_lda( rs, rx, ry, rz, v0, v1, f0, f1, f2, 1.0,
                                       ts, tx, ty, tz, ks, kx, ky, kz );
      else
        xckernel::nc_fxc_contract_limit_lda( rs, v0, v1, f0, f1, f2,
                                             ts, tx, ty, tz, ks, kx, ky, kz );
      const double kslot[4] = { ks, kz, ky, kx };
      for( int k = 0; k < 4; ++k ) {
        const double c = 0.5 * weights[ip] * kslot[k];
        const auto* b = basis_eval.data() + size_t(ip)*nbe;
        auto* z = Z.data() + k*blk + size_t(ip)*nbe;
        for( int32_t mu = 0; mu < nbe; ++mu ) z[mu] = c * b[mu];
      }
    }

    for( int k = 0; k < 4; ++k )
      lwd->inc_vxc( npts, nbf, nbe, basis_eval.data(), submat_map, Z.data() + k*blk, nbe,
                    FXC[k], ldFXC[k], nbe_scr.data() );
  } // tasks
  } // omp parallel

  for( int k = 0; k < 4; ++k )
    for( int32_t j = 0; j < nbf; ++j )
      for( int32_t i = j+1; i < nbf; ++i )
        FXC[k][ j + i*ldFXC[k] ] = FXC[k][ i + j*ldFXC[k] ];
}

} // namespace GauXC::detail
