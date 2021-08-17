#include "local_host_work_driver_pimpl.hpp"

namespace GauXC {

struct ReferenceLocalHostWorkDriver : public detail::LocalHostWorkDriverPIMPL {

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

  void eval_xmat( size_t npts, size_t nbf, size_t nbe, 
    const submat_map_t& submat_map, const double* P, size_t ldp, 
    const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) 
    override;
    
  void eval_uvvar_lda( size_t npts, size_t nbe, const double* basis_eval,
    const double* X, size_t ldx, double* den_eval) override;
  void eval_uvvar_gga( size_t npts, size_t nbe, const double* basis_eval,
    const double* dbasis_x_eval, const double *dbasis_y_eval, 
    const double* dbasis_z_eval, const double* X, size_t ldx, double* den_eval, 
    double* dden_x_eval, double* dden_y_eval, double* dden_z_eval, 
    double* gamma ) override;

  void eval_zmat_lda_vxc( size_t npts, size_t nbe, const double* vrho, 
    const double* basis_eval, double* Z, size_t ldz ) override;
  void eval_zmat_gga_vxc( size_t npts, size_t nbe, const double* vrho, 
    const double* vgamma, const double* basis_eval, const double* dbasis_x_eval,
    const double* dbasis_y_eval, const double* dbasis_z_eval, 
    const double* dden_x_eval, const double* dden_y_eval, const double* dden_z_eval,
    double* Z, size_t ldz ) override;


  void inc_vxc( size_t npts, size_t nbf, size_t nbe, 
    const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
    size_t ldz, double* VXC, size_t ldvxc, double* scr ) override;

};

}