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

  virtual void eval_xmat( size_t npts, size_t nbf, size_t nbe, 
    const submat_map_t& submat_map, const double* P, size_t ldp, 
    const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) = 0;
    
  virtual void eval_uvvar_lda( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) = 0;
  virtual void eval_uvvar_gga( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eavl, const double *dbasis_y_eval, 
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval, 
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval, 
    double* gamma ) = 0;

  virtual void eval_zmat_lda_vxc( size_t npts, size_t nbe, const double* vrho, 
    const double* basis_eval, double* Z, size_t ldz ) = 0;
  virtual void eval_zmat_gga_vxc( size_t npts, size_t nbe, const double* vrho, 
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