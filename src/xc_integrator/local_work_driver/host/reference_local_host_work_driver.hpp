/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include "local_host_work_driver_pimpl.hpp"

namespace GauXC {

struct ReferenceLocalHostWorkDriver : public detail::LocalHostWorkDriverPIMPL {

  double *boys_table;
  
  using submat_map_t   = LocalHostWorkDriverPIMPL::submat_map_t;
  using task_container = LocalHostWorkDriverPIMPL::task_container;
  using tast_iterator  = LocalHostWorkDriverPIMPL::task_iterator;

  ReferenceLocalHostWorkDriver();

  virtual ~ReferenceLocalHostWorkDriver() noexcept;

  ReferenceLocalHostWorkDriver( const ReferenceLocalHostWorkDriver& )     = delete;
  ReferenceLocalHostWorkDriver( ReferenceLocalHostWorkDriver&& ) noexcept = delete;

  // Public APIs

  void partition_weights( XCWeightAlg weight_alg, const Molecule& mol, 
    const MolMeta& meta, task_iterator task_begin, task_iterator task_end ) override;

  void eval_collocation( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval ) override;
  void eval_collocation_gradient( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval, double* dbasis_x_eval, double* dbasis_y_eval, 
    double* dbasis_z_eval) override;
  void eval_collocation_hessian( size_t npts, size_t nshells, size_t nbe, 
    const double* pts, const BasisSet<double>& basis, const int32_t* shell_list, 
    double* basis_eval, double* dbasis_x_eval, double* dbasis_y_eval, 
    double* dbasis_z_eval, double* d2basis_xx_eval, double* d2basis_xy_eval,
    double* d2basis_xz_eval, double* d2basis_yy_eval, double* d2basis_yz_eval,
    double* d2basis_zz_eval ) override;

  void eval_xmat( size_t npts, size_t nbf, size_t nbe, 
    const submat_map_t& submat_map, const double* P, size_t ldp, 
    const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) 
    override;

  void eval_exx_gmat( size_t npts, size_t nshells, size_t nshell_pairs,
    size_t nbe, const double* points, const double* weights, 
    const BasisSet<double>& basis, const ShellPairCollection<double>& shpairs, 
    const BasisSetMap& basis_map, const int32_t* shell_list, 
    const std::pair<int32_t,int32_t>* shell_pair_list, 
    const double* X, size_t ldx, double* G, size_t ldg ) override ;

  void eval_exx_fmat( size_t npts, size_t nbf, size_t nbe_bra,
    size_t nbe_ket, const submat_map_t& submat_map_bra,
    const submat_map_t& submat_map_ket, const double* P, size_t ldp,
    const double* basis_eval, size_t ldb, double* F, size_t ldf,
    double* scr ) override;

  void inc_exx_k( size_t npts, size_t nbf, size_t nbe_bra, size_t nbe_ket, 
    const double* basis_eval, const submat_map_t& submat_map_bra, 
    const submat_map_t& submat_map_ket, const double* G, size_t ldg, double* K, 
    size_t ldk, double* scr ) override;
    
  void eval_uvvar_lda_rks( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) override;

  void eval_uvvar_lda_uks( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) override;
  void eval_uvvar_gga_rks( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eval, const double *dbasis_y_eval, 
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval, 
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval, 
    double* gamma ) override;
  void eval_uvvar_gga_uks( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eval, const double *dbasis_y_eval,
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval,
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval,
    double* gamma ) override;


  void eval_zmat_lda_vxc_rks( size_t npts, size_t nbe, const double* vrho, 
    const double* basis_eval, double* Z, size_t ldz ) override;
  void eval_zmat_lda_vxc_uks( size_t npts, size_t nbe, const double* vrho,
    const double* basis_eval, double* Z, size_t ldz ) override;


  void eval_zmat_gga_vxc_rks( size_t npts, size_t nbe, const double* vrho, 
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval, 
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Z, size_t ldz ) override;
  void eval_zmat_gga_vxc_uks( size_t npts, size_t nbe, const double* vrho,
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval,
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Z, size_t ldz ) override;




  void inc_vxc( size_t npts, size_t nbf, size_t nbe, 
    const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
    size_t ldz, double* VXC, size_t ldvxc, double* scr ) override;

};

}
