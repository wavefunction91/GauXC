#include "host/reference_local_host_work_driver.hpp"
#include "host/reference/weights.hpp"
#include "host/reference/collocation.hpp"

#include "host/util.hpp"
#include "host/blas.hpp"
#include <stdexcept>

#include <gauxc/basisset_map.hpp>
#include "rys_integral.h"
#include <gauxc/util/real_solid_harmonics.hpp>

namespace GauXC {

ReferenceLocalHostWorkDriver::ReferenceLocalHostWorkDriver() = default; 
ReferenceLocalHostWorkDriver::~ReferenceLocalHostWorkDriver() noexcept = default;

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
      throw std::runtime_error("Weight Alg Not Supported");
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
  const submat_map_t& submat_map, const double* P, size_t ldp, 
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

  blas::gemm( 'N', 'N', nbe, npts, nbe, 2., P_use, ldp_use, basis_eval, ldb, 
    0., X, ldx );

}




// U/VVar LDA (density)
void ReferenceLocalHostWorkDriver::eval_uvvar_lda( size_t npts, size_t nbe, 
 const double* basis_eval, const double* X, size_t ldx, double* den_eval) {

  for( int32_t i = 0; i < (int32_t)npts; ++i ) {

    const size_t ioff = size_t(i) * ldx;
    const auto*   X_i = X + ioff;
    den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );

  }

}



// U/VVar GGA (density + grad, gamma)
void ReferenceLocalHostWorkDriver::eval_uvvar_gga( size_t npts, size_t nbe, 
  const double* basis_eval, const double* dbasis_x_eval, 
  const double *dbasis_y_eval, const double* dbasis_z_eval, const double* X, 
  size_t ldx, double* den_eval, double* dden_x_eval, double* dden_y_eval, 
  double* dden_z_eval, double* gamma ) {

  for( int32_t i = 0; i < (int32_t)npts; ++i ) {

    const size_t ioff = size_t(i) * ldx;
    const auto*   X_i = X + ioff;

    den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );

    const auto dx = 
      2. * blas::dot( nbe, dbasis_x_eval + ioff, 1, X_i, 1 );
    const auto dy = 
      2. * blas::dot( nbe, dbasis_y_eval + ioff, 1, X_i, 1 );
    const auto dz = 
      2. * blas::dot( nbe, dbasis_z_eval + ioff, 1, X_i, 1 );

    dden_x_eval[i] = dx;
    dden_y_eval[i] = dy;
    dden_z_eval[i] = dz;

    gamma[i] = dx*dx + dy*dy + dz*dz;

  }

}







// Eval Z Matrix LDA VXC
void ReferenceLocalHostWorkDriver::eval_zmat_lda_vxc( size_t npts, size_t nbf, 
  const double* vrho, const double* basis_eval, double* Z, size_t ldz ) {


  blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Z, ldz );

  for( int32_t i = 0; i < (int32_t)npts; ++i ) {

    auto* z_col = Z + i*ldz;

    const double fact = 0.5 * vrho[i];
    GauXC::blas::scal( nbf, fact, z_col, 1 );

  }

}



// Eval Z Matrix GGA VXC
void ReferenceLocalHostWorkDriver::eval_zmat_gga_vxc( size_t npts, size_t nbf, 
  const double* vrho, const double* vgamma, const double* basis_eval, 
  const double* dbasis_x_eval, const double* dbasis_y_eval, 
  const double* dbasis_z_eval, const double* dden_x_eval, 
  const double* dden_y_eval, const double* dden_z_eval, double* Z, size_t ldz ) {

  if( ldz != nbf ) 
    throw std::logic_error(std::string("INVALID DIMS")+std::string(__PRETTY_FUNCTION__));
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






struct RysBasis {
  std::vector< shells > _shells;
  RysBasis( const BasisSet<double>& basis ) {
    size_t nshells = basis.size();
    _shells.resize(nshells);
    for( int i = 0; i < nshells; ++i ) {
      _shells[i].origin.x = basis[i].O()[0];
      _shells[i].origin.y = basis[i].O()[1];
      _shells[i].origin.z = basis[i].O()[2];

      _shells[i].m = basis[i].nprim();
      _shells[i].L = basis[i].l();

      _shells[i].coeff = new coefficients[_shells[i].m];
      for( int j = 0; j < _shells[i].m; ++j ) {
        _shells[i].coeff[j].alpha = basis[i].alpha()[j];
        _shells[i].coeff[j].coeff = basis[i].coeff()[j];
      }
    }
  }

  shells& operator[](int i){ return _shells.at(i); }
  const shells& operator[](int i) const { return _shells.at(i); }

  ~RysBasis() noexcept {
    for( auto& sh : _shells ) delete sh.coeff;
  }
};






void ReferenceLocalHostWorkDriver::
  eval_exx_gmat( size_t npts, size_t nbe, const double* points, const double* weights, 
    const BasisSet<double>& basis, const BasisSetMap& basis_map, const double* X, size_t ldx, 
    double* G, size_t ldg ) {

  // Cast points to Rys format (binary compatable)
  point* _points = reinterpret_cast<point*>(const_cast<double*>(points));

  const size_t nshells = basis.nshells();

  // Set G to zero
  for( int j = 0; j < npts; ++j )
  for( int i = 0; i < nbe;  ++i ) {
    G[i + j*ldg] = 0.;
  }

  // Copy the basis set 
  RysBasis rys_basis(basis);

  // Spherical Harmonic Transformer
  util::SphericalHarmonicTransform sph_trans(5);

  size_t ioff = 0;
  for( int i = 0; i < nshells; ++i ) {
    const auto& bra       = basis.at(i);
    const int bra_l       = bra.l();
    const int bra_sz      = bra.size();
    const int bra_cart_sz = (bra_l+1) * (bra_l+2) / 2;
    const bool need_transform_bra = bra.pure() and bra_l > 0;

    size_t joff = 0;
    for( int j = 0; j <= i; ++j ) {
      const auto& ket       = basis.at(j);
      const int ket_l       = ket.l();
      const int ket_sz      = ket.size();
      const int ket_cart_sz = (ket_l+1) * (ket_l+2) / 2;
      const bool need_transform_ket = ket.pure() and ket_l > 0;

      //const int shpair_sz      = bra_sz * ket_sz;
      const int shpair_cart_sz = bra_cart_sz * ket_cart_sz;

      std::vector<double> _tmp( shpair_cart_sz * npts );

      compute_integral_shell_pair( npts, rys_basis[i], rys_basis[j], _points,
        _tmp.data() );

      const auto need_tform = need_transform_bra or need_transform_ket;

      const auto* Xj = X + joff;
      const auto* Xi = X + ioff;
      auto*       Gj = G + joff;
      auto*       Gi = G + ioff;
      for( int k = 0; k < npts; ++k ) {
        auto Xjk = Xj + k*ldx;
        auto Xik = Xi + k*ldx;
        auto Gjk = Gj + k*ldg;
        auto Gik = Gi + k*ldg;

        std::vector<double> _tmp_sph( bra_sz * ket_sz, 0. );
        auto* ints = _tmp.data() + k * shpair_cart_sz;
        if( need_transform_bra and need_transform_ket )
          sph_trans.tform_both( bra_l, ket_l, ints, ket_cart_sz, 
            _tmp_sph.data(), ket_sz );
        else if( need_transform_bra )
          sph_trans.tform_bra( bra_l, ket_sz, ints, ket_cart_sz, 
            _tmp_sph.data(), ket_sz );
        else if( need_transform_ket )
          sph_trans.tform_ket( bra_sz, ket_l, ints, ket_cart_sz, 
            _tmp_sph.data(), ket_sz );

        if( need_tform ) ints = _tmp_sph.data();

        for( int ii = 0; ii < bra_sz; ++ii )
        for( int jj = 0; jj < ket_sz; ++jj ) {
          Gik[ii] += weights[k] * ints[ii*ket_sz + jj] * Xjk[jj];
	      }

       	if( i != j )
        for( int ii = 0; ii < bra_sz; ++ii )
        for( int jj = 0; jj < ket_sz; ++jj ) {
          Gjk[jj] += weights[k] * ints[ii*ket_sz + jj] * Xik[ii];
	      }
      }

      joff += ket_sz;
    }
    ioff += bra_sz;
  }


}








}
