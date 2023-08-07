/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "local_host_work_driver.hpp"

namespace GauXC  {
namespace detail {

struct LocalHostWorkDriverPIMPL {

  using submat_map_t   = LocalHostWorkDriver::submat_map_t;
  using task_container = LocalHostWorkDriver::task_container;
  using task_iterator  = LocalHostWorkDriver::task_iterator;

  LocalHostWorkDriverPIMPL();

  virtual ~LocalHostWorkDriverPIMPL() noexcept;

  LocalHostWorkDriverPIMPL( const LocalHostWorkDriverPIMPL& )     = delete;
  LocalHostWorkDriverPIMPL( LocalHostWorkDriverPIMPL&& ) noexcept = delete;


  // Public APIs

  virtual void partition_weights( XCWeightAlg weight_alg, const Molecule& mol, 
    const MolMeta& meta, task_iterator task_begin, task_iterator task_end ) = 0;

  virtual void eval_collocation( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval ) = 0;
  virtual void eval_collocation_gradient( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval, double* dbasis_x_eval, double* dbasis_y_eval, 
    double* dbasis_z_eval) = 0;
  virtual void eval_collocation_hessian( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval, double* dbasis_x_eval, double* dbasis_y_eval, 
    double* dbasis_z_eval, double* d2basis_xx_eval, double* d2basis_xy_eval,
    double* d2basis_xz_eval, double* d2basis_yy_eval, double* d2basis_yz_eval,
    double* d2basis_zz_eval ) = 0;

  virtual void eval_xmat( size_t npts, size_t nbf, size_t nbe, 
    const submat_map_t& submat_map, const double* P, size_t ldp, 
    const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) = 0;

  virtual void eval_exx_fmat( size_t npts, size_t nbf, size_t nbe_bra,
    size_t nbe_ket, const submat_map_t& submat_map_bra,
    const submat_map_t& submat_map_ket, const double* P, size_t ldp,
    const double* basis_eval, size_t ldb, double* F, size_t ldf,
    double* scr ) = 0;

  virtual void eval_exx_gmat( size_t npts, size_t nshells, size_t nshell_pairs,
    size_t nbe, const double* points, const double* weights, 
    const BasisSet<double>& basis, const ShellPairCollection<double>& shpairs, 
    const BasisSetMap& basis_map, const int32_t* shell_list, 
    const std::pair<int32_t,int32_t>* shell_pair_list, 
    const double* X, size_t ldx, double* G, size_t ldg ) = 0;

  virtual void inc_exx_k( size_t npts, size_t nbf, size_t nbe_bra, size_t nbe_ket, 
    const double* basis_eval, const submat_map_t& submat_map_bra, 
    const submat_map_t& submat_map_ket, const double* G, size_t ldg, double* K, 
    size_t ldk, double* scr ) = 0;
    
  virtual void eval_uvvar_lda_rks( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) = 0;
  virtual void eval_uvvar_lda_uks( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) = 0;

  virtual void eval_uvvar_gga_rks( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eavl, const double *dbasis_y_eval, 
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval, 
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval, 
    double* gamma ) = 0;
  virtual void eval_uvvar_gga_uks( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eavl, const double *dbasis_y_eval,
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval,
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval,
    double* gamma ) = 0;


  virtual void eval_zmat_lda_vxc_rks( size_t npts, size_t nbe, const double* vrho, 
    const double* basis_eval, double* Z, size_t ldz ) = 0;
  virtual void eval_zmat_lda_vxc_uks( size_t npts, size_t nbe, const double* vrho,
    const double* basis_eval, double* Z, size_t ldz ) = 0;

  virtual void eval_zmat_gga_vxc_rks( size_t npts, size_t nbe, const double* vrho, 
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval, 
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Z, size_t ldz ) = 0;
  virtual void eval_zmat_gga_vxc_uks( size_t npts, size_t nbe, const double* vrho,
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval,
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Z, size_t ldz ) = 0;

  virtual void inc_vxc( size_t npts, size_t nbf, size_t nbe, 
    const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
    size_t ldz, double* VXC, size_t ldvxc, double* scr ) = 0;

};


}
}
