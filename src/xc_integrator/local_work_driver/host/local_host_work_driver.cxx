/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "local_host_work_driver_pimpl.hpp"
#include <stdexcept>

namespace GauXC {

LocalHostWorkDriver::LocalHostWorkDriver() : 
  pimpl_(nullptr) { }
LocalHostWorkDriver::LocalHostWorkDriver(pimpl_type&& ptr) :
  pimpl_( std::move(ptr) ){ }

LocalHostWorkDriver::~LocalHostWorkDriver() noexcept = default;

LocalHostWorkDriver::LocalHostWorkDriver( LocalHostWorkDriver&& other ) noexcept :
  pimpl_(std::move(other.pimpl_)) { }

#define throw_if_invalid_pimpl(ptr) \
  if(not ptr) GAUXC_PIMPL_NOT_INITIALIZED()





// Partition weights
void LocalHostWorkDriver::partition_weights( XCWeightAlg weight_alg, 
  const Molecule& mol, const MolMeta& meta, task_iterator task_begin, 
  task_iterator task_end ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->partition_weights(weight_alg, mol, meta, task_begin, task_end);

}


// Collocation
void LocalHostWorkDriver::eval_collocation( size_t npts, size_t nshells, size_t nbe, 
  const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
  double* basis_eval ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_collocation(npts, nshells, nbe, pts, basis, shell_list, basis_eval);

}


// Collocation Gradient
void LocalHostWorkDriver::eval_collocation_gradient( size_t npts, size_t nshells, 
  size_t nbe, const double* pts, const BasisSet<double>& basis, 
  const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
  double* dbasis_y_eval, double* dbasis_z_eval) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_collocation_gradient(npts, nshells, nbe, pts, basis, shell_list, basis_eval,
    dbasis_x_eval, dbasis_y_eval, dbasis_z_eval);

}


// Collocation Hessian
void LocalHostWorkDriver::eval_collocation_hessian( size_t npts, size_t nshells, 
    size_t nbe, const double* pts, const BasisSet<double>& basis, 
    const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
    double* dbasis_y_eval, double* dbasis_z_eval, double* d2basis_xx_eval, 
    double* d2basis_xy_eval, double* d2basis_xz_eval, double* d2basis_yy_eval, 
    double* d2basis_yz_eval, double* d2basis_zz_eval ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_collocation_hessian(npts, nshells, nbe, pts, basis, shell_list, basis_eval,
    dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval, d2basis_xy_eval,
    d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval, d2basis_zz_eval);

}


// X matrix (P * B)
void LocalHostWorkDriver::eval_xmat( size_t npts, size_t nbf, size_t nbe, 
  const submat_map_t& submat_map, const double* P, size_t ldp, 
  const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_xmat(npts, nbf, nbe, submat_map, P, ldp, basis_eval, ldb, X, 
    ldx, scr);

}


void LocalHostWorkDriver::eval_exx_fmat( size_t npts, size_t nbf, size_t nbe_bra,
  size_t nbe_ket, const submat_map_t& submat_map_bra,
  const submat_map_t& submat_map_ket, const double* P, size_t ldp,
  const double* basis_eval, size_t ldb, double* F, size_t ldf,
  double* scr ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_exx_fmat(npts, nbf, nbe_bra, nbe_ket, submat_map_bra,
    submat_map_ket, P, ldp, basis_eval, ldb, F, ldf, scr ); 

}


// G Matrix G(mu,i) = w(i) * A(mu,nu,i) * X(mu,i)
void LocalHostWorkDriver::eval_exx_gmat( size_t npts, size_t nshells, 
  size_t nshell_pairs, size_t nbe, const double* points, const double* weights, 
  const BasisSet<double>& basis, const ShellPairCollection<double>& shpairs, 
  const BasisSetMap& basis_map, const int32_t* shell_list, 
  const std::pair<int32_t,int32_t>* shell_pair_list, 
  const double* X, size_t ldx, double* G, size_t ldg ) {;

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_exx_gmat(npts, nshells, nshell_pairs, nbe, points, weights, 
    basis, shpairs, basis_map, shell_list, shell_pair_list, X, ldx, G, ldg );

}

void LocalHostWorkDriver::inc_exx_k( size_t npts, size_t nbf, size_t nbe_bra, 
  size_t nbe_ket, const double* basis_eval, const submat_map_t& submat_map_bra, 
  const submat_map_t& submat_map_ket, const double* G, size_t ldg, double* K, 
  size_t ldk, double* scr ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->inc_exx_k(npts, nbf, nbe_bra, nbe_ket, basis_eval, submat_map_bra,
    submat_map_ket, G, ldg, K, ldk, scr );
}



// U/VVar LDA (density)
void LocalHostWorkDriver::eval_uvvar_lda( size_t npts, size_t nbe, 
 const double* basis_eval, const double* X, size_t ldx, double* den_eval) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_uvvar_lda(npts, nbe, basis_eval, X, ldx, den_eval);

}


// U/VVar GGA (density + grad, gamma)
void LocalHostWorkDriver::eval_uvvar_gga( size_t npts, size_t nbe, 
  const double* basis_eval, const double* dbasis_x_eval, 
  const double *dbasis_y_eval, const double* dbasis_z_eval, const double* X, 
  size_t ldx, double* den_eval, double* dden_x_eval, double* dden_y_eval, 
  double* dden_z_eval, double* gamma ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_uvvar_gga(npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
    dbasis_z_eval, X, ldx, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
    gamma);

}

// Eval Z Matrix LDA VXC
void LocalHostWorkDriver::eval_zmat_lda_vxc( size_t npts, size_t nbe, 
  const double* vrho, const double* basis_eval, double* Z, size_t ldz ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_zmat_lda_vxc(npts, nbe, vrho, basis_eval, Z, ldz);

}

// Eval Z Matrix GGA VXC
void LocalHostWorkDriver::eval_zmat_gga_vxc( size_t npts, size_t nbe, 
  const double* vrho, const double* vgamma, const double* basis_eval, 
  const double* dbasis_x_eval, const double* dbasis_y_eval, 
  const double* dbasis_z_eval, const double* dden_x_eval, 
  const double* dden_y_eval, const double* dden_z_eval, double* Z, size_t ldz ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->eval_zmat_gga_vxc(npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
    dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval, dden_z_eval,
    Z, ldz);

}


// Increment VXC by Z
void LocalHostWorkDriver::inc_vxc( size_t npts, size_t nbf, size_t nbe, 
  const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
  size_t ldz, double* VXC, size_t ldvxc, double* scr ) {

  throw_if_invalid_pimpl(pimpl_);
  pimpl_->inc_vxc(npts, nbf, nbe, basis_eval, submat_map, Z, ldz, VXC, ldvxc,
    scr);

}



}
