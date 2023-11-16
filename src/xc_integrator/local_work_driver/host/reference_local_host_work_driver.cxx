/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */

#include "host/reference_local_host_work_driver.hpp"
#include "host/reference/weights.hpp"
#include "host/reference/collocation.hpp"

#include "host/util.hpp"
#include "host/blas.hpp"
#include <stdexcept>

#include <gauxc/basisset_map.hpp>
#include <gauxc/shell_pair.hpp>
#include <gauxc/util/unused.hpp>
#include "cpu/integral_data_types.hpp"
#include "cpu/obara_saika_integrals.hpp"
#include "cpu/chebyshev_boys_computation.hpp"
#include <gauxc/util/real_solid_harmonics.hpp>
#include "integrator_util/integral_bounds.hpp"

namespace GauXC {

  ReferenceLocalHostWorkDriver::ReferenceLocalHostWorkDriver() {
    this->boys_table = XCPU::boys_init();
  }
  
  ReferenceLocalHostWorkDriver::~ReferenceLocalHostWorkDriver() noexcept {
    XCPU::boys_finalize(this->boys_table);
  }

  // Partition weights
  void ReferenceLocalHostWorkDriver::partition_weights( XCWeightAlg weight_alg, 
							const Molecule& mol, const MolMeta& meta, task_iterator task_begin, 
							task_iterator task_end ) {
    switch( weight_alg ) {
      case XCWeightAlg::Becke:
        reference_becke_weights_host( mol, meta, task_begin, task_end );
        break;
      case XCWeightAlg::SSF:
        reference_ssf_weights_host( mol, meta, task_begin, task_end );
        break;
      case XCWeightAlg::LKO:
        reference_lko_weights_host( mol, meta, task_begin, task_end );
        break;
      default:
        GAUXC_GENERIC_EXCEPTION("Weight Alg Not Supported");
    }
  }


  // Collocation
  void ReferenceLocalHostWorkDriver::eval_collocation( size_t npts, size_t nshells, 
						       size_t nbe, const double* pts, const BasisSet<double>& basis, 
						       const int32_t* shell_list, double* basis_eval ) {
    gau2grid_collocation( npts, nshells, nbe, pts, basis, shell_list, basis_eval );
  }


  // Collocation Gradient
  void ReferenceLocalHostWorkDriver::eval_collocation_gradient( size_t npts, 
								size_t nshells, size_t nbe, const double* pts, const BasisSet<double>& basis, 
								const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
								double* dbasis_y_eval, double* dbasis_z_eval) {
    gau2grid_collocation_gradient(npts, nshells, nbe, pts, basis, shell_list,
				  basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
  }

  void ReferenceLocalHostWorkDriver::eval_collocation_hessian( size_t npts, 
							       size_t nshells, size_t nbe, const double* pts, const BasisSet<double>& basis, 
							       const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
							       double* dbasis_y_eval, double* dbasis_z_eval, double* d2basis_xx_eval, 
							       double* d2basis_xy_eval, double* d2basis_xz_eval, double* d2basis_yy_eval, 
							       double* d2basis_yz_eval, double* d2basis_zz_eval ) {
    gau2grid_collocation_hessian(npts, nshells, nbe, pts, basis, shell_list,
				 basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval,
				 d2basis_xy_eval, d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval,
				 d2basis_zz_eval);
  }


  // X matrix (P * B)
  void ReferenceLocalHostWorkDriver::eval_xmat( size_t npts, size_t nbf, size_t nbe, 
						const submat_map_t& submat_map, double fac, const double* P, size_t ldp, 
						const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) {
    const auto* P_use = P;
    size_t ldp_use = ldp;
     
    if( submat_map.size() > 1 ) {
      detail::submat_set( nbf, nbf, nbe, nbe, P, ldp, scr, nbe, submat_map );
      P_use = scr;
      ldp_use = nbe;
    } else if( nbe != nbf ) {
      P_use = P + submat_map[0][0]*(ldp+1);
    }

    blas::gemm( 'N', 'N', nbe, npts, nbe, fac, P_use, ldp_use, basis_eval, ldb, 
		0., X, ldx );

  }

  // U/VVar LDA (density)
  void ReferenceLocalHostWorkDriver::eval_uvvar_lda_rks( size_t npts, size_t nbe, 
						     const double* basis_eval, const double* X, size_t ldx, double* den_eval) {


    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioff = size_t(i) * ldx;
      const auto*   X_i = X + ioff;
      den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );

    }    

  }

  
  void ReferenceLocalHostWorkDriver::eval_uvvar_lda_uks( size_t npts, size_t nbe,
   const double* basis_eval, const double* Xs, size_t ldxs, 
   const double* Xz, size_t ldxz, double* den_eval) {
  
    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioffs = size_t(i) * ldxs;
      const size_t ioffz = size_t(i) * ldxz;

      const auto*   Xs_i = Xs + ioffs;
      const auto*   Xz_i = Xz + ioffz;

      const double rhos = blas::dot( nbe, basis_eval + ioffs, 1, Xs_i, 1 );
      const double rhoz = blas::dot( nbe, basis_eval + ioffz, 1, Xz_i, 1 );
      
      den_eval[2*i]   = 0.5*(rhos + rhoz); // rho_+
      den_eval[2*i+1] = 0.5*(rhos - rhoz); // rho_-

    }
 
  }
  
  void ReferenceLocalHostWorkDriver::eval_uvvar_lda_gks( size_t npts, size_t nbe, const double* basis_eval,
    const double* Xs, size_t ldxs, const double* Xz, size_t ldxz,
    const double* Xx, size_t ldxx, const double* Xy, size_t ldxy, double* den_eval, double* K) {


    auto *K2 = K; // KZ // store K in the Z matrix
    auto *K3 = K2 + npts;
    auto *K4 = K3 + npts;
 
    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioffs = size_t(i) * ldxs;
      const size_t ioffz = size_t(i) * ldxz;
      const size_t ioffx = size_t(i) * ldxx;
      const size_t ioffy = size_t(i) * ldxy;

      const auto*   Xs_i = Xs + ioffs;
      const auto*   Xz_i = Xz + ioffz;
      const auto*   Xx_i = Xx + ioffx;
      const auto*   Xy_i = Xy + ioffy;

      const double rhos = blas::dot( nbe, basis_eval + ioffs, 1, Xs_i, 1 );
      const double rhoz = blas::dot( nbe, basis_eval + ioffz, 1, Xz_i, 1 );
      const double rhox = blas::dot( nbe, basis_eval + ioffx, 1, Xx_i, 1 );
      const double rhoy = blas::dot( nbe, basis_eval + ioffy, 1, Xy_i, 1 );
 
      double mtemp = rhoz * rhoz + rhox * rhox + rhoy * rhoy;
      double mnorm = 0;

      if (mtemp > 1.0e-24) {
        mnorm = sqrt(mtemp);
        K2[i] = rhoz / mnorm;
        K3[i] = rhox / mnorm;
        K4[i] = rhoy / mnorm;
      } else {
        mnorm = (1. / 3.) * (rhox + rhoy + rhoz);
        K2[i] = 1. / 3.;
        K3[i] = 1. / 3.;
        K4[i] = 1. / 3.;
      }

      den_eval[2*i]   = 0.5*(rhos + mnorm); // rho_+
      den_eval[2*i+1] = 0.5*(rhos - mnorm); // rho_-

    }

  }


  void ReferenceLocalHostWorkDriver::eval_uvvar_gga_rks( size_t npts, size_t nbe, 
						     const double* basis_eval, const double* dbasis_x_eval, 
						     const double *dbasis_y_eval, const double* dbasis_z_eval, const double* X, 
						     size_t ldx, double* den_eval, double* dden_x_eval, double* dden_y_eval, 
						     double* dden_z_eval, double* gamma ) {

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioff = size_t(i) * ldx;
      const auto*   X_i = X + ioff;

      den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );

      const auto dx = 2. * blas::dot( nbe, dbasis_x_eval + ioff, 1, X_i, 1 );
      const auto dy = 2. * blas::dot( nbe, dbasis_y_eval + ioff, 1, X_i, 1 );
      const auto dz = 2. * blas::dot( nbe, dbasis_z_eval + ioff, 1, X_i, 1 );

      dden_x_eval[i] = dx;
      dden_y_eval[i] = dy;
      dden_z_eval[i] = dz;

      gamma[i] = dx*dx + dy*dy + dz*dz;

    }
  }

void ReferenceLocalHostWorkDriver::eval_uvvar_gga_uks( size_t npts, size_t nbe,
  const double* basis_eval, const double* dbasis_x_eval,
  const double *dbasis_y_eval, const double* dbasis_z_eval, const double* Xs,
  size_t ldxs, const double* Xz, size_t ldxz, 
  double* den_eval, double* dden_x_eval, double* dden_y_eval,
  double* dden_z_eval, double* gamma ) {

   for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioffs = size_t(i) * ldxs;
      const size_t ioffz = size_t(i) * ldxz;

      const auto*   Xs_i = Xs + ioffs;
      const auto*   Xz_i = Xz + ioffz;

      double rhos = blas::dot( nbe, basis_eval + ioffs, 1, Xs_i, 1 ); // S density
      double rhoz = blas::dot( nbe, basis_eval + ioffz, 1, Xz_i, 1 ); // Z density


      den_eval[2*i]   = 0.5*(rhos + rhoz); // rho_+
      den_eval[2*i+1] = 0.5*(rhos - rhoz); // rho_-

      const auto dndx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffs, 1, Xs_i, 1 );
      const auto dndy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffs, 1, Xs_i, 1 );
      const auto dndz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffs, 1, Xs_i, 1 );

      const auto dMzdx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffz, 1, Xz_i, 1 );
      const auto dMzdy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffz, 1, Xz_i, 1 );
      const auto dMzdz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffz, 1, Xz_i, 1 );

      dden_x_eval[2*i] = dndx; // dn / dx
      dden_y_eval[2*i] = dndy; // dn / dy
      dden_z_eval[2*i] = dndz; // dn / dz

      dden_x_eval[2*i+1] = dMzdx; // dMz / dx
      dden_y_eval[2*i+1] = dMzdy; // dMz / dy
      dden_z_eval[2*i+1] = dMzdz; // dMz / dz

      // (del n).(del n)
      const auto dn_sq  = dndx*dndx + dndy*dndy + dndz*dndz;
      // (del Mz).(del Mz)
      const auto dMz_sq = dMzdx*dMzdx + dMzdy*dMzdy + dMzdz*dMzdz;
      // (del n).(del Mz)
      const auto dn_dMz = dndx*dMzdx + dndy*dMzdy + dndz*dMzdz;

      gamma[3*i  ] = 0.25*(dn_sq + dMz_sq) + 0.5*dn_dMz;
      gamma[3*i+1] = 0.25*(dn_sq - dMz_sq);
      gamma[3*i+2] = 0.25*(dn_sq + dMz_sq) - 0.5*dn_dMz;
    }

}

void ReferenceLocalHostWorkDriver::eval_uvvar_gga_gks( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eval, const double *dbasis_y_eval,
    const double* dbasis_z_eval, const double* Xs, size_t ldxs,
    const double* Xz, size_t ldxz, const double* Xx, size_t ldxx,
    const double* Xy, size_t ldxy, double* den_eval,
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval, double* gamma, double* K, double* H ) {

   auto *K2 = K; // KZ // store K in the Z matrix
   auto *K3 = K2 + npts;
   auto *K4 = K3 + npts;

   auto *H2 = H; // KZ // store K in the Z matrix
   auto *H3 = H2 + npts;
   auto *H4 = H3 + npts;

   for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioffs = size_t(i) * ldxs;
      const size_t ioffz = size_t(i) * ldxz;
      const size_t ioffx = size_t(i) * ldxx;
      const size_t ioffy = size_t(i) * ldxy;

      const auto*   Xs_i = Xs + ioffs;
      const auto*   Xz_i = Xz + ioffz;
      const auto*   Xx_i = Xx + ioffx;
      const auto*   Xy_i = Xy + ioffy;

      const double rhos = blas::dot( nbe, basis_eval + ioffs, 1, Xs_i, 1 );
      const double rhoz = blas::dot( nbe, basis_eval + ioffz, 1, Xz_i, 1 );
      const double rhox = blas::dot( nbe, basis_eval + ioffx, 1, Xx_i, 1 );
      const double rhoy = blas::dot( nbe, basis_eval + ioffy, 1, Xy_i, 1 );

      const auto dndx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffs, 1, Xs_i, 1 );
      const auto dndy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffs, 1, Xs_i, 1 );
      const auto dndz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffs, 1, Xs_i, 1 );

      const auto dMzdx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffz, 1, Xz_i, 1 );
      const auto dMzdy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffz, 1, Xz_i, 1 );
      const auto dMzdz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffz, 1, Xz_i, 1 );

      const auto dMxdx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffx, 1, Xx_i, 1 );
      const auto dMxdy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffx, 1, Xx_i, 1 );
      const auto dMxdz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffx, 1, Xx_i, 1 );

      const auto dMydx =
        2. * blas::dot( nbe, dbasis_x_eval + ioffy, 1, Xy_i, 1 );
      const auto dMydy =
        2. * blas::dot( nbe, dbasis_y_eval + ioffy, 1, Xy_i, 1 );
      const auto dMydz =
        2. * blas::dot( nbe, dbasis_z_eval + ioffy, 1, Xy_i, 1 );


      dden_x_eval[4 * i] = dndx;
      dden_y_eval[4 * i] = dndy;
      dden_z_eval[4 * i] = dndz;

      dden_x_eval[4 * i + 1] = dMzdx;
      dden_y_eval[4 * i + 1] = dMzdy;
      dden_z_eval[4 * i + 1] = dMzdz;

      dden_x_eval[4 * i + 2] = dMxdx;
      dden_y_eval[4 * i + 2] = dMxdy;
      dden_z_eval[4 * i + 2] = dMxdz;

      dden_x_eval[4 * i + 3] = dMydx;
      dden_y_eval[4 * i + 3] = dMydy;
      dden_z_eval[4 * i + 3] = dMydz;

      double mtemp = rhoz * rhoz + rhox * rhox + rhoy * rhoy;
      double mnorm = 0;

      auto dels_dot_dels = dndx * dndx + dndy * dndy + dndz * dndz;
      auto delz_dot_delz = dMzdx * dMzdx + dMzdy * dMzdy + dMzdz * dMzdz;
      auto delx_dot_delx = dMxdx * dMxdx + dMxdy * dMxdy + dMxdz * dMxdz;
      auto dely_dot_dely = dMydx * dMydx + dMydy * dMydy + dMydz * dMydz;

      auto dels_dot_delz = dndx * dMzdx + dndy * dMzdy + dndz * dMzdz;
      auto dels_dot_delx = dndx * dMxdx + dndy * dMxdy + dndz * dMxdz;
      auto dels_dot_dely = dndx * dMydx + dndy * dMydy + dndz * dMydz;

      auto sum = delz_dot_delz + delx_dot_delx + dely_dot_dely;
      auto s_sum =
          dels_dot_delz * rhoz + dels_dot_delx * rhox + dels_dot_dely * rhoy;

      auto sqsum2 =
          sqrt(dels_dot_delz * dels_dot_delz + dels_dot_delx * dels_dot_delx +
               dels_dot_dely * dels_dot_dely);

      double sign = 1.;
      if (std::signbit(s_sum))
        sign = -1.;

      if (mtemp > 1.0e-24) {
        mnorm = sqrt(mtemp);
        K2[i] = rhoz / mnorm;
        K3[i] = rhox / mnorm;
        K4[i] = rhoy / mnorm;
        H2[i] = sign * dels_dot_delz / sqsum2;
        H3[i] = sign * dels_dot_delx / sqsum2;
        H4[i] = sign * dels_dot_dely / sqsum2;
      } else {
        mnorm = (1. / 3.) * (rhox + rhoy + rhoz);
        K2[i] = 1. / 3.;
        K3[i] = 1. / 3.;
        K4[i] = 1. / 3.;

        H2[i] = sign / 3.;
        H3[i] = sign / 3.;
        H4[i] = sign / 3.;
      }
      
      den_eval[2 * i] = 0.5 * (rhos + mnorm);
      den_eval[2 * i + 1] = 0.5 * (rhos - mnorm);
      
      gamma[3 * i] = 0.25 * (dels_dot_dels + sum) + 0.5 * sign * sqsum2;
      gamma[3 * i + 1] = 0.25 * (dels_dot_dels - sum);
      gamma[3 * i + 2] = 0.25 * (dels_dot_dels + sum) - 0.5 * sign * sqsum2;


    }

}

  // Eval Z Matrix LDA VXC
  void ReferenceLocalHostWorkDriver::eval_zmat_lda_vxc_rks( size_t npts, size_t nbf, 
							const double* vrho, const double* basis_eval, double* Z, size_t ldz ) {


    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Z, ldz );

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      auto* z_col = Z + i*ldz;

      const double fact = 0.5 * vrho[i];
      GauXC::blas::scal( nbf, fact, z_col, 1 );

    }

  }

  // Eval Z Matrix LDA VXC
  void ReferenceLocalHostWorkDriver::eval_zmat_lda_vxc_uks( size_t npts, size_t nbf,
              const double* vrho, const double* basis_eval, double* Zs, size_t ldzs,
              double* Zz, size_t ldzz ) {


    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zs, ldzs);
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zz, ldzz);

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      auto* zs_col = Zs + i*ldzs;
      auto* zz_col = Zz + i*ldzz;

      const double factp = 0.5 * vrho[2*i];
      const double factm = 0.5 * vrho[2*i+1];

      //eq. 56 https://doi.org/10.1140/epjb/e2018-90170-1
      GauXC::blas::scal( nbf, 0.5*(factp + factm), zs_col, 1 );
      GauXC::blas::scal( nbf, 0.5*(factp - factm), zz_col, 1 );

    }
 

  }

void ReferenceLocalHostWorkDriver::eval_zmat_lda_vxc_gks( size_t npts, size_t nbe, const double* vrho,
    const double* basis_eval, double* Zs, size_t ldzs, double* Zz, size_t ldzz,
    double* Zx, size_t ldzx,double* Zy, size_t ldzy, double *K ) {

  auto *KZ = K; // KZ // store K in the Z matrix
  auto *KX = KZ + npts;
  auto *KY = KX + npts;

    blas::lacpy( 'A', nbe, npts, basis_eval, nbe, Zs, ldzs);
    blas::lacpy( 'A', nbe, npts, basis_eval, nbe, Zz, ldzz);
    blas::lacpy( 'A', nbe, npts, basis_eval, nbe, Zx, ldzx);
    blas::lacpy( 'A', nbe, npts, basis_eval, nbe, Zy, ldzy);

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      auto* zs_col = Zs + i*ldzs;
      auto* zz_col = Zz + i*ldzz;
      auto* zx_col = Zx + i*ldzx;
      auto* zy_col = Zy + i*ldzy;

      const double factp = 0.5 * vrho[2*i];
      const double factm = 0.5 * vrho[2*i+1];
      const double factor = 0.5 * (factp - factm);

      //eq. 56 https://doi.org/10.1140/epjb/e2018-90170-1
      GauXC::blas::scal( nbe, 0.5*(factp + factm), zs_col, 1 );
      GauXC::blas::scal( nbe, KZ[i] * factor, zz_col, 1 );
      GauXC::blas::scal( nbe, KX[i] * factor, zx_col, 1 );
      GauXC::blas::scal( nbe, KY[i] * factor, zy_col, 1 );
   
    }

}

  // Eval Z Matrix GGA VXC
  void ReferenceLocalHostWorkDriver::eval_zmat_gga_vxc_rks( size_t npts, size_t nbf, 
							const double* vrho, const double* vgamma, const double* basis_eval, 
							const double* dbasis_x_eval, const double* dbasis_y_eval, 
							const double* dbasis_z_eval, const double* dden_x_eval, 
							const double* dden_y_eval, const double* dden_z_eval, double* Z, size_t ldz ) {

    if( ldz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Z, nbf );

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const int32_t ioff = i * nbf;

      auto* z_col    = Z + ioff;
      auto* bf_x_col = dbasis_x_eval + ioff; 
      auto* bf_y_col = dbasis_y_eval + ioff; 
      auto* bf_z_col = dbasis_z_eval + ioff; 

      const auto lda_fact = 0.5 * vrho[i];
      blas::scal( nbf, lda_fact, z_col, 1 );

      const auto gga_fact = 2. * vgamma[i]; 
      const auto x_fact = gga_fact * dden_x_eval[i];
      const auto y_fact = gga_fact * dden_y_eval[i];
      const auto z_fact = gga_fact * dden_z_eval[i];

      blas::axpy( nbf, x_fact, bf_x_col, 1, z_col, 1 );
      blas::axpy( nbf, y_fact, bf_y_col, 1, z_col, 1 );
      blas::axpy( nbf, z_fact, bf_z_col, 1, z_col, 1 );

    }

  }

  void ReferenceLocalHostWorkDriver::eval_zmat_gga_vxc_uks( size_t npts, size_t nbf,
              const double* vrho, const double* vgamma, const double* basis_eval,
              const double* dbasis_x_eval, const double* dbasis_y_eval,
              const double* dbasis_z_eval, const double* dden_x_eval,
              const double* dden_y_eval, const double* dden_z_eval, double* Zs, 
              size_t ldzs, double* Zz, size_t ldzz ) {


    if( ldzs != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    if( ldzz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zs, ldzs);
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zz, ldzz);

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const int32_t ioff = i * nbf;

      auto* zs_col = Zs + ioff;
      auto* zz_col = Zz + ioff;
      auto* bf_x_col = dbasis_x_eval + ioff;
      auto* bf_y_col = dbasis_y_eval + ioff;
      auto* bf_z_col = dbasis_z_eval + ioff;

      const double factp = 0.5 * vrho[2*i];
      const double factm = 0.5 * vrho[2*i+1];

      GauXC::blas::scal( nbf, 0.5*(factp + factm), zs_col, 1 ); //additional 0.5 is from eq 56 in petrone 2018 eur phys journal b "an efficent implementation of .. "
      GauXC::blas::scal( nbf, 0.5*(factp - factm), zz_col, 1 );

      const auto gga_fact_pp = vgamma[3*i];
      const auto gga_fact_pm = vgamma[3*i+1];
      const auto gga_fact_mm = vgamma[3*i+2];

      const auto gga_fact_1 = 0.5*(gga_fact_pp + gga_fact_pm + gga_fact_mm);
      const auto gga_fact_2 = 0.5*(gga_fact_pp - gga_fact_mm);
      const auto gga_fact_3 = 0.5*(gga_fact_pp - gga_fact_pm + gga_fact_mm);

      const auto x_fact_s = gga_fact_1 * dden_x_eval[2*i] + gga_fact_2 * dden_x_eval[2*i+1];
      const auto y_fact_s = gga_fact_1 * dden_y_eval[2*i] + gga_fact_2 * dden_y_eval[2*i+1];
      const auto z_fact_s = gga_fact_1 * dden_z_eval[2*i] + gga_fact_2 * dden_z_eval[2*i+1];

      const auto x_fact_z = gga_fact_3 * dden_x_eval[2*i+1] + gga_fact_2 * dden_x_eval[2*i];
      const auto y_fact_z = gga_fact_3 * dden_y_eval[2*i+1] + gga_fact_2 * dden_y_eval[2*i];
      const auto z_fact_z = gga_fact_3 * dden_z_eval[2*i+1] + gga_fact_2 * dden_z_eval[2*i];
      
      blas::axpy( nbf, x_fact_s, bf_x_col, 1, zs_col, 1 );
      blas::axpy( nbf, y_fact_s, bf_y_col, 1, zs_col, 1 );
      blas::axpy( nbf, z_fact_s, bf_z_col, 1, zs_col, 1 );

      blas::axpy( nbf, x_fact_z, bf_x_col, 1, zz_col, 1 );
      blas::axpy( nbf, y_fact_z, bf_y_col, 1, zz_col, 1 );
      blas::axpy( nbf, z_fact_z, bf_z_col, 1, zz_col, 1 );

    }
  }

void ReferenceLocalHostWorkDriver::eval_zmat_gga_vxc_gks( size_t npts, size_t nbf, const double* vrho,
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval,
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Zs, size_t ldzs, double* Zz, size_t ldzz, double* Zx, size_t ldzx,
    double* Zy, size_t ldzy, double* K, double* H ) {

    auto *KZ = K; // KZ // store K in the Z matrix
    auto *KX = KZ + npts;
    auto *KY = KX + npts;

    auto *HZ = H; // KZ // store K in the Z matrix
    auto *HX = HZ + npts;
    auto *HY = HX + npts;

    if( ldzs != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    if( ldzz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    if( ldzx != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    if( ldzy != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));

    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zs, ldzs);
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zz, ldzz);
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zx, ldzx);
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Zy, ldzy);   


    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const int32_t ioff = i * nbf;

      auto* zs_col = Zs + ioff;
      auto* zz_col = Zz + ioff;
      auto* zx_col = Zx + ioff;
      auto* zy_col = Zy + ioff;

      auto* bf_x_col = dbasis_x_eval + ioff;
      auto* bf_y_col = dbasis_y_eval + ioff;
      auto* bf_z_col = dbasis_z_eval + ioff;

      const double factp = 0.5 * vrho[2*i];
      const double factm = 0.5 * vrho[2*i+1];
      const double factor = 0.5 * (factp - factm);

      GauXC::blas::scal( nbf, 0.5*(factp + factm), zs_col, 1 ); //additional 0.5 is from eq 56 in petrone 2018 eur phys journal b "an efficent implementation of .. "
      GauXC::blas::scal( nbf, KZ[i]*factor, zz_col, 1 );
      GauXC::blas::scal( nbf, KX[i]*factor, zx_col, 1 );
      GauXC::blas::scal( nbf, KY[i]*factor, zy_col, 1 );

      const auto gga_fact_pp = vgamma[3 * i];
      const auto gga_fact_pm = vgamma[3 * i + 1];
      const auto gga_fact_mm = vgamma[3 * i + 2];

      const auto gga_fact_1 = 0.5 * (gga_fact_pp + gga_fact_pm + gga_fact_mm);
      const auto gga_fact_2 = 0.5 * (gga_fact_pp - gga_fact_mm);
      const auto gga_fact_3 = 0.5 * (gga_fact_pp - gga_fact_pm + gga_fact_mm);

      const auto x_fact_s = gga_fact_1 * dden_x_eval[4 * i] +
                            gga_fact_2 * (HZ[i] * dden_x_eval[4 * i + 1] +
                                          HX[i] * dden_x_eval[4 * i + 2] +
                                          HY[i] * dden_x_eval[4 * i + 3]);
      const auto y_fact_s = gga_fact_1 * dden_y_eval[4 * i] +
                            gga_fact_2 * (HZ[i] * dden_y_eval[4 * i + 1] +
                                          HX[i] * dden_y_eval[4 * i + 2] +
                                          HY[i] * dden_y_eval[4 * i + 3]);
      const auto z_fact_s = gga_fact_1 * dden_z_eval[4 * i] +
                            gga_fact_2 * (HZ[i] * dden_z_eval[4 * i + 1] +
                                          HX[i] * dden_z_eval[4 * i + 2] +
                                          HY[i] * dden_z_eval[4 * i + 3]);

      const auto x_fact_z = gga_fact_3 * dden_x_eval[4 * i + 1] +
                            gga_fact_2 * HZ[i] * dden_x_eval[4 * i];
      const auto y_fact_z = gga_fact_3 * dden_y_eval[4 * i + 1] +
                            gga_fact_2 * HZ[i] * dden_y_eval[4 * i];
      const auto z_fact_z = gga_fact_3 * dden_z_eval[4 * i + 1] +
                            gga_fact_2 * HZ[i] * dden_z_eval[4 * i];

      const auto x_fact_x = gga_fact_3 * dden_x_eval[4 * i + 2] +
                            gga_fact_2 * HX[i] * dden_x_eval[4 * i];
      const auto y_fact_x = gga_fact_3 * dden_y_eval[4 * i + 2] +
                            gga_fact_2 * HX[i] * dden_y_eval[4 * i];
      const auto z_fact_x = gga_fact_3 * dden_z_eval[4 * i + 2] +
                            gga_fact_2 * HX[i] * dden_z_eval[4 * i];

      const auto x_fact_y = gga_fact_3 * dden_x_eval[4 * i + 3] +
                            gga_fact_2 * HY[i] * dden_x_eval[4 * i];
      const auto y_fact_y = gga_fact_3 * dden_y_eval[4 * i + 3] +
                            gga_fact_2 * HY[i] * dden_y_eval[4 * i];
      const auto z_fact_y = gga_fact_3 * dden_z_eval[4 * i + 3] +
                            gga_fact_2 * HY[i] * dden_z_eval[4 * i];


      blas::axpy(nbf, x_fact_s, bf_x_col, 1, zs_col, 1);
      blas::axpy(nbf, y_fact_s, bf_y_col, 1, zs_col, 1);
      blas::axpy(nbf, z_fact_s, bf_z_col, 1, zs_col, 1);

      blas::axpy(nbf, x_fact_z, bf_x_col, 1, zz_col, 1);
      blas::axpy(nbf, y_fact_z, bf_y_col, 1, zz_col, 1);
      blas::axpy(nbf, z_fact_z, bf_z_col, 1, zz_col, 1);

      blas::axpy(nbf, x_fact_x, bf_x_col, 1, zx_col, 1);
      blas::axpy(nbf, y_fact_x, bf_y_col, 1, zx_col, 1);
      blas::axpy(nbf, z_fact_x, bf_z_col, 1, zx_col, 1);

      blas::axpy(nbf, x_fact_y, bf_x_col, 1, zy_col, 1);
      blas::axpy(nbf, y_fact_y, bf_y_col, 1, zy_col, 1);
      blas::axpy(nbf, z_fact_y, bf_z_col, 1, zy_col, 1);

    }

}

  // Increment VXC by Z
  void ReferenceLocalHostWorkDriver::inc_vxc( size_t npts, size_t nbf, size_t nbe, 
					      const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
					      size_t ldz, double* VXC, size_t ldvxc, double* scr ) {

    if( submat_map.size() > 1 ) {
      blas::syr2k('L', 'N', nbe, npts, 1., basis_eval, nbe, Z, ldz, 0., scr, nbe );
      detail::inc_by_submat( nbf, nbf, nbe, nbe, VXC, ldvxc, scr, nbe, submat_map );
    } else {
      blas::syr2k('L', 'N', nbe, npts, 1., basis_eval, nbe, Z, ldz, 1., 
		  VXC + submat_map[0][0]*(ldvxc+1), ldvxc );
    }

  }

  // Increment K by G
  void ReferenceLocalHostWorkDriver::inc_exx_k( size_t npts, size_t nbf, 
						size_t nbe_bra, size_t nbe_ket, const double* basis_eval, 
						const submat_map_t& submat_map_bra, const submat_map_t& submat_map_ket, 
						const double* G, size_t ldg, double* K, size_t ldk, double* scr ) {

    if( submat_map_bra.size() > 1 or submat_map_ket.size() > 1 ) {
      blas::gemm( 'N', 'T', nbe_bra, nbe_ket, npts, 1., basis_eval, nbe_bra,
		  G, ldg, 0., scr, nbe_bra );

      detail::inc_by_submat( nbf, nbf, nbe_bra, nbe_ket, K, ldk, scr, nbe_bra, 
			     submat_map_bra, submat_map_ket );
    } else {
      blas::gemm( 'N', 'T', nbe_bra, nbe_ket, npts, 1., basis_eval, nbe_bra,
		  G, ldg, 1., K + submat_map_ket[0][0]*ldk + submat_map_bra[0][0], ldk );
    }

  }


  // Construct F = P * B (P non-square, TODO: should merge with XMAT)
  void ReferenceLocalHostWorkDriver::eval_exx_fmat( size_t npts, size_t nbf, 
						    size_t nbe_bra, size_t nbe_ket, const submat_map_t& submat_map_bra,
						    const submat_map_t& submat_map_ket, const double* P, size_t ldp,
						    const double* basis_eval, size_t ldb, double* F, size_t ldf,
						    double* scr ) {

    const auto* P_use = P;
    size_t ldp_use = ldp;

    if( submat_map_bra.size() > 1 or submat_map_ket.size() > 1 ) {
      detail::submat_set( nbf, nbf, nbe_bra, nbe_ket, P, ldp,
			  scr, nbe_bra, submat_map_bra, submat_map_ket );
      P_use = scr;
      ldp_use = nbe_bra;
    } else {
      P_use = P + submat_map_ket[0][0]*ldp + submat_map_bra[0][0];
    }

    blas::gemm( 'N', 'N', nbe_bra, npts, nbe_ket, 1., P_use, ldp_use, basis_eval,
		ldb, 0., F, ldf );

  }

  // Construct G(mu,i) = w(i) * A(mu,nu,i) * F(nu, i)
  void ReferenceLocalHostWorkDriver::eval_exx_gmat( size_t npts, size_t nshells, 
    size_t nshell_pairs, size_t nbe, const double* points, const double* weights, 
    const BasisSet<double>& basis, const ShellPairCollection<double>& shpairs, 
    const BasisSetMap& basis_map, const int32_t* shell_list, 
    const std::pair<int32_t,int32_t>* shell_pair_list, 
    const double* X, size_t ldx, double* G, size_t ldg ) {

    util::unused(basis_map);

    // Cast points to Rys format (binary compatable)
    XCPU::point* _points = 
      reinterpret_cast<XCPU::point*>(const_cast<double*>(points));
    std::vector<double> _points_transposed(3 * npts);

    for(size_t i = 0; i < npts; ++i) {
      _points_transposed[i + 0 * npts] = _points[i].x;
      _points_transposed[i + 1 * npts] = _points[i].y;
      _points_transposed[i + 2 * npts] = _points[i].z;
    }

  
    // Set G to zero
    for( size_t j = 0; j < npts; ++j )
    for( size_t i = 0; i < nbe;  ++i ) {
	    G[i + j*ldg] = 0.;
    }


    // Spherical Harmonic Transformer
    util::SphericalHarmonicTransform sph_trans(5);

    const bool any_pure = std::any_of( shell_list, shell_list + nshells,
				       [&](const auto& i){ return basis.at(i).pure(); } );
    
    const size_t nbe_cart = 
      basis.nbf_cart_subset( shell_list, shell_list + nshells );

    std::vector<double> X_cart, G_cart;
    if( any_pure ){
      X_cart.resize( nbe_cart * npts );
      G_cart.resize( nbe_cart * npts, 0. );

      // Transform X into cartesian
      int ioff = 0;
      int ioff_cart = 0;
      for( auto i = 0ul; i < nshells; ++i ) {
        const auto ish = shell_list[i];
        const auto& shell      = basis.at(ish);
        const int shell_l       = shell.l();
        const int shell_sz      = shell.size();
        const int shell_cart_sz = shell.cart_size();
        
        if( shell.pure() and shell_l > 0 ) {
          sph_trans.itform_bra_cm( shell_l, npts, X + ioff, ldx,
        			   X_cart.data() + ioff_cart, nbe_cart );
        } else {
          blas::lacpy( 'A', shell_sz, npts, X + ioff, ldx,
        	       X_cart.data() + ioff_cart, nbe_cart );
        }
        ioff += shell_sz;
        ioff_cart += shell_cart_sz;
      }
    }

    const auto* X_use = any_pure ? X_cart.data() : X;
    auto*       G_use = any_pure ? G_cart.data() : G;
    const auto ldx_use = any_pure ? nbe_cart : ldx;
    const auto ldg_use = any_pure ? nbe_cart : ldg;

    std::vector<double> X_cart_rm( nbe_cart*npts,0. ), 
                        G_cart_rm( nbe_cart*npts,0. );
    for( auto i = 0ul; i < nbe_cart; ++i )
    for( auto j = 0ul; j < npts;     ++j ) {
      X_cart_rm[i*npts + j] = X_use[i + j*ldx_use];
    }


    std::map<size_t,size_t> cou_offsets_map;
    std::vector<size_t> cou_cart_sizes(nshells);
    cou_cart_sizes[0] = 0;
    cou_offsets_map[shell_list[0]] = 0;
    for(size_t i = 1; i < nshells; ++i) {
      cou_cart_sizes[i] = cou_cart_sizes[i-1] +
        basis.at(shell_list[i-1]).cart_size();
      cou_offsets_map[shell_list[i]] = cou_cart_sizes[i];
    }

    size_t ndo = 0;
    {
#if 0
    //size_t ioff_cart = 0;
    for( auto i = 0ul; i < nshells; ++i ) {
      const auto ish        = shell_list[i];
      const auto& bra       = basis[ish];
      const int bra_cart_sz = bra.cart_size();
      const size_t ioff_cart = cou_cart_sizes[i] * npts;
      XCPU::point bra_origin{bra.O()[0],bra.O()[1],bra.O()[2]};

      //size_t joff_cart = 0;
      for( auto j = 0ul; j <= i; ++j ) {
      //for( auto j = i; j < nshells; ++j ) {
        const auto jsh        = shell_list[j];
        const auto& ket       = basis[jsh];
        const int ket_cart_sz = ket.cart_size();
        const size_t joff_cart = cou_cart_sizes[j] * npts;
        XCPU::point ket_origin{ket.O()[0],ket.O()[1],ket.O()[2]};
        if(!need_sp(ish,jsh)) continue;
        ++ndo;

        auto sh_pair = shpairs.at(ish,jsh);
        auto prim_pair_data = sh_pair.prim_pairs();
        auto nprim_pair     = sh_pair.nprim_pairs();
        
        XCPU::compute_integral_shell_pair( ish == jsh,
        				   npts, _points_transposed.data(),
        				   bra.l(), ket.l(), bra_origin, ket_origin,
        				   nprim_pair, prim_pair_data,
        				   X_cart_rm.data()+ioff_cart, X_cart_rm.data()+joff_cart, npts,
        				   G_cart_rm.data()+ioff_cart, G_cart_rm.data()+joff_cart, npts,
        				   const_cast<double*>(weights), this->boys_table );
        
        //joff_cart += ket_cart_sz * npts;
      }
	
      //ioff_cart += bra_cart_sz * npts;
    }
#else
    for( auto ij = 0ul; ij < nshell_pairs; ++ij ) {
      auto [ish,jsh] = shell_pair_list[ij];
      //std::cout << "SHP " << ij << " " << i << " " << j << " " << nshells << std::endl;

     
      // Bra
      const auto& bra      = basis.at(ish);
      const auto ioff_cart = cou_offsets_map.at(ish) * npts;
      XCPU::point bra_origin{bra.O()[0],bra.O()[1],bra.O()[2]};

      // Ket
      const auto& ket      = basis.at(jsh);
      const auto joff_cart = cou_offsets_map.at(jsh) * npts;
      XCPU::point ket_origin{ket.O()[0],ket.O()[1],ket.O()[2]};

      auto sh_pair = shpairs.at(ish,jsh);
      auto prim_pair_data = sh_pair.prim_pairs();
      auto nprim_pair     = sh_pair.nprim_pairs();
      
      ndo++;  
      XCPU::compute_integral_shell_pair( ish == jsh,
      				   npts, _points_transposed.data(),
      				   bra.l(), ket.l(), bra_origin, ket_origin,
      				   nprim_pair, prim_pair_data,
      				   X_cart_rm.data()+ioff_cart, X_cart_rm.data()+joff_cart, npts,
      				   G_cart_rm.data()+ioff_cart, G_cart_rm.data()+joff_cart, npts,
      				   const_cast<double*>(weights), this->boys_table );
    }
#endif
    }
    //std::cout << "NDO " << ndo << " " << ndo / double(nshells*(nshells+1)/2) << std::endl;
   
    for( auto i = 0ul; i < nbe_cart; ++i )
    for( auto j = 0ul; j < npts;     ++j ) {
	    G_use[i + j*ldg_use] = G_cart_rm[i*npts + j];
    }
  
    // Transform G back to spherical
    if( any_pure ) {
      size_t ioff = 0;
      size_t ioff_cart = 0;
      for( auto i = 0ul; i < nshells; ++i ) {
        const auto ish = shell_list[i];
        const auto& shell      = basis.at(ish);
        const int shell_l       = shell.l();
        const int shell_sz      = shell.size();
        const int shell_cart_sz = shell.cart_size();
        
        if( shell.pure() and shell_l > 0 ) {
          sph_trans.tform_bra_cm( shell_l, npts, G_cart.data() + ioff_cart, nbe_cart,
        			  G + ioff, ldg );
        } else {
          blas::lacpy( 'A', shell_sz, npts, G_cart.data() + ioff_cart, nbe_cart,
        	       G + ioff, ldg );
        }
        ioff += shell_sz;
        ioff_cart += shell_cart_sz;
      }
    }

  } // GMAT

}
