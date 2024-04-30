/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_task.hpp>
#include <gauxc/util/div_ceil.hpp>
#include <vector>
#include <gauxc/basisset_map.hpp>
#include <gauxc/shell_pair.hpp>
#include <gauxc/molmeta.hpp>
//#include <gauxc/reduction_driver.hpp>
#include <any>
#include <cstring>
#include "device/device_queue.hpp"

namespace GauXC {

enum integrator_xc_approx : uint32_t {
  _UNDEFINED = 0,
  LDA        = 1,
  GGA        = 2,
  MGGA_TAU   = 3,
  MGGA_LAPL  = 4
};

struct integrator_term_tracker {
  bool weights                   = false;
  bool den                       = false;
  bool exc_vxc                   = false;
  bool exc_grad                  = false;
  bool exx                       = false;
  bool exx_ek_screening          = false;
  integrator_xc_approx xc_approx = _UNDEFINED;
  inline void reset() {
    std::memset( this, 0, sizeof(integrator_term_tracker) );
  }
};

#define PRDVL(pred,val) (pred) ? (val) : 0ul;

struct required_term_storage {
  bool grid_points  = false;
  bool grid_weights = false;

  inline size_t grid_points_size(size_t npts) { 
    return PRDVL(grid_points, 3 * npts); 
  }
  inline size_t grid_weights_size(size_t npts) { 
    return PRDVL(grid_weights, npts); 
  }

  // Evaluation of functions on the grid (linear storage)
  bool grid_den      = false;
  bool grid_den_grad = false;
  bool grid_den_lapl = false;
  bool grid_gamma    = false;
  bool grid_tau      = false;
  bool grid_eps      = false;
  bool grid_vrho     = false;
  bool grid_vgamma   = false;
  bool grid_vtau     = false;
  bool grid_vlapl    = false;

  inline size_t grid_den_size(size_t npts){ 
    return PRDVL(grid_den, npts);
  }
  inline size_t grid_den_grad_size(size_t npts){ 
    return PRDVL(grid_den_grad, 3 * npts);
  }
  inline size_t grid_den_lapl_size(size_t npts){ 
    return PRDVL(grid_den_lapl, npts);
  }
  inline size_t grid_gamma_size(size_t npts){ 
    return PRDVL(grid_gamma, npts);
  }
  inline size_t grid_tau_size(size_t npts){ 
    return PRDVL(grid_tau, npts);
  }
  inline size_t grid_eps_size(size_t npts){ 
    return PRDVL(grid_eps, npts);
  }
  inline size_t grid_vrho_size(size_t npts){ 
    return PRDVL(grid_vrho, npts);
  }
  inline size_t grid_vgamma_size(size_t npts){ 
    return PRDVL(grid_vgamma, npts);
  }
  inline size_t grid_vtau_size(size_t npts){ 
    return PRDVL(grid_vtau, npts);
  }
  inline size_t grid_vlapl_size(size_t npts){ 
    return PRDVL(grid_vlapl, npts);
  }

  // Task-local matrices
  bool task_bfn           = false;
  bool task_bfn_grad      = false;
  bool task_bfn_hess      = false;
  bool task_bfn_lapl      = false;
  bool task_zmat          = false;
  bool task_xmat          = false;
  bool task_xmat_grad     = false;
  bool task_fmat          = false;
  bool task_gmat          = false;
  bool task_nbe_scr       = false;
  bool task_bfn_shell_indirection = false;


  inline size_t task_bfn_size(size_t nbe, size_t npts) {
    return PRDVL(task_bfn, nbe * npts);
  }
  inline size_t task_bfn_grad_size(size_t nbe, size_t npts) {
    return PRDVL(task_bfn_grad, 3 * nbe * npts);
  }
  inline size_t task_bfn_hess_size(size_t nbe, size_t npts) {
    return PRDVL(task_bfn_hess, 6 * nbe * npts);
  }
  inline size_t task_bfn_lapl_size(size_t nbe, size_t npts) {
    return PRDVL(task_bfn_lapl, nbe * npts);
  }
  inline size_t task_zmat_size(size_t nbe, size_t npts) {
    return PRDVL(task_zmat, nbe * npts);
  }
  inline size_t task_xmat_grad_size(size_t nbe, size_t npts) {
    return PRDVL(task_xmat_grad, 3 * nbe * npts);
  }
  inline size_t task_fmat_size(size_t nbe, size_t npts) {
    return PRDVL(task_fmat, nbe * npts);
  }
  inline size_t task_gmat_size(size_t nbe, size_t npts) {
    return PRDVL(task_gmat, nbe * npts);
  }
  inline size_t task_nbe_scr_size(size_t nbe_bfn, size_t nbe_cou) {
    return PRDVL(task_nbe_scr, std::max(nbe_bfn,nbe_cou) * nbe_bfn);
  }
  inline size_t task_bfn_shell_indirection_size(size_t nbe) {
    return PRDVL(task_bfn_shell_indirection, nbe);
  }

  // Index packing
  bool task_submat_cut_bfn   = false;
  bool task_submat_block_bfn = false;
  bool task_submat_cut_cou   = false;
  bool task_submat_block_cou = false;

  inline size_t task_submat_cut_bfn_size(size_t ncut) {
    return PRDVL(task_submat_cut_bfn, 3*ncut);
  }
  inline size_t task_submat_block_bfn_size(size_t nblock) {
    return PRDVL(task_submat_block_bfn, nblock);
  }
  inline size_t task_submat_cut_cou_size(size_t ncut) {
    return PRDVL(task_submat_cut_cou, 3*ncut);
  }
  inline size_t task_submat_block_cou_size(size_t nblock) {
    return PRDVL(task_submat_block_cou, nblock);
  }

  // Task indirection
  bool task_indirection = false;
  inline size_t task_indirection_size() {
    return PRDVL(task_indirection, 1ul);
  }

  // Weights kernel scratch
  bool grid_to_center_dist_scr     = false;
  bool grid_to_center_dist_nearest = false;
  bool grid_to_parent_center       = false;

  inline size_t grid_to_center_dist_scr_size(size_t ldatom, size_t npts) {
    return PRDVL(grid_to_center_dist_scr, ldatom * npts);
  }
  inline size_t grid_to_center_dist_nearest_size(size_t npts) {
    return PRDVL(grid_to_center_dist_nearest, npts);
  }
  inline size_t grid_to_parent_center_size(size_t npts) {
    return PRDVL(grid_to_parent_center, npts);
  }

  // Shell/Shell pairs lists + indirection
  bool task_shell_list_bfn    = false;
  bool task_shell_offs_bfn    = false;
  bool shell_to_task_bfn      = false;
  bool shell_pair_to_task_cou = false;
  bool task_to_shell_pair_cou = false;
  
  inline size_t task_shell_list_bfn_size(size_t nshells) {
    return PRDVL(task_shell_list_bfn, nshells);
  }
  inline size_t task_shell_offs_bfn_size(size_t nshells) {
    return PRDVL(task_shell_offs_bfn, nshells);
  }
  inline size_t shell_to_task_idx_bfn_size(size_t nshells) {
    return PRDVL(shell_to_task_bfn, nshells);
  }
  inline size_t shell_to_task_off_bfn_size(size_t nshells) {
    return PRDVL(shell_to_task_bfn, nshells);
  }
  inline size_t shell_pair_to_task_idx_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(shell_pair_to_task_cou, nslt);
  }
  inline size_t shell_pair_to_task_row_off_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(shell_pair_to_task_cou, nslt);
  }
  inline size_t shell_pair_to_task_col_off_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(shell_pair_to_task_cou, nslt);
  }
  inline size_t task_to_shell_pair_col_off_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(task_to_shell_pair_cou, nslt);
  }
  inline size_t task_to_shell_pair_row_off_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(task_to_shell_pair_cou, nslt);
  }
  inline size_t task_to_shell_pair_idx_cou_size(size_t nshells) {
    const size_t nslt = (nshells * (nshells+1)) / 2;
    return PRDVL(task_to_shell_pair_cou, nslt);
  }
  inline size_t task_to_shell_pair_cou_size() {
    return PRDVL(task_to_shell_pair_cou, 1ul);
  }
  inline size_t task_to_shell_pair_cou_subtask_size(size_t npts, size_t subtask_size) {
    const size_t num_subtasks = util::div_ceil(npts, subtask_size);
    return PRDVL(task_to_shell_pair_cou, num_subtasks);
  }



  inline explicit required_term_storage(integrator_term_tracker tracker) {
    // Everything under the sun needs the grid
    grid_points  = true;
    grid_weights = true;

    if(tracker.weights) {
      grid_to_center_dist_scr     = true;
      grid_to_center_dist_nearest = true; 
      grid_to_parent_center       = true;
    }

    // Allocated terms for XC aclculations
    const bool is_xc = tracker.exc_vxc or tracker.exc_grad;
    if(is_xc) {
      if( tracker.xc_approx == _UNDEFINED )
        GAUXC_GENERIC_EXCEPTION("NO XC APPROX SET");
      //const bool is_lda  = is_xc and tracker.xc_approx == LDA;
      const bool is_gga  = is_xc and tracker.xc_approx == GGA;
      const bool need_tau  = tracker.xc_approx == MGGA_TAU;
      const bool need_lapl = tracker.xc_approx == MGGA_LAPL;
      const bool is_mgga = is_xc and (need_tau or need_lapl);
      const bool is_grad = tracker.exc_grad;

      grid_den      = true;
      grid_den_grad = is_gga or is_mgga or is_grad;
      grid_den_lapl = need_lapl;
      grid_gamma    = is_gga or is_mgga;
      grid_tau      = is_mgga;
      grid_eps      = true;
      grid_vrho     = true;
      grid_vgamma   = is_gga or is_mgga;
      grid_vtau     = is_mgga;
      grid_vlapl    = need_lapl;

      task_bfn          = true;
      task_bfn_grad     = is_gga or  is_mgga or is_grad;
      task_bfn_hess     = is_gga and is_grad;
      task_bfn_lapl     = need_lapl;
      task_zmat         = true;
      task_xmat         = true;
      task_xmat_grad    = is_mgga or (is_gga and is_grad);
      task_nbe_scr      = true;

      task_submat_cut_bfn   = true;
      task_submat_block_bfn = true;
      task_indirection      = true;

      task_shell_list_bfn = true;
      task_shell_offs_bfn = true;
      shell_to_task_bfn   = true;
    }

    // Density integration
    if(tracker.den) {
      grid_den              = true;
      task_bfn              = true;
      task_nbe_scr          = true;
      task_xmat             = true;
      task_zmat             = true;
      task_submat_cut_bfn   = true;
      task_submat_block_bfn = true;
      task_indirection      = true;

      task_shell_list_bfn = true;
      task_shell_offs_bfn = true;
      shell_to_task_bfn   = true;
    }

    // EXX integration
    if(tracker.exx) {
      task_bfn              = true;
      task_fmat             = true;
      task_gmat             = true;
      task_nbe_scr          = true;
      task_submat_cut_bfn   = true;
      task_submat_block_bfn = true;
      task_submat_cut_cou   = true;
      task_submat_block_cou = true;
      task_indirection      = true;

      task_shell_list_bfn    = true;
      task_shell_offs_bfn    = true;
      shell_to_task_bfn      = true;
      //shell_pair_to_task_cou = true;
      task_to_shell_pair_cou = true;
    }

    if(tracker.exx_ek_screening) {
      task_bfn              = true;
      task_indirection      = true;

      task_shell_list_bfn        = true;
      task_shell_offs_bfn        = true;
      task_bfn_shell_indirection = true;
      shell_to_task_bfn          = true;
    }

  }
};

#undef PRDVL



inline 
std::ostream& operator<<( std::ostream& out, const integrator_term_tracker& t ) {
  out << std::boolalpha;
  out << "Integrator Terms:" << std::endl;
  out << "  WEIGHTS  " << t.weights << std::endl;
  out << "  DEN      " << t.den << std::endl;
  out << "  EXC_VXC  " << t.exc_vxc << std::endl;
  out << "  EXC_GRAD " << t.exc_grad << std::endl;
  out << "  EXX      " << t.exx << std::endl;
  return out;
}

/** Base class for all XCDeviceData types
 *
 *  Exposes virtual API to manage device memory and batch XC
 *  integration tasks.
 */
struct XCDeviceData {

  using host_task_type        = XCTask;
  using host_task_container   = std::vector<host_task_type>;
  using host_task_iterator    = host_task_container::iterator;

  virtual ~XCDeviceData() noexcept = default;

  /// Allocate device memory for data that will persist on the device.
  virtual void reset_allocations() = 0;
  virtual void allocate_static_data_weights( int32_t natoms ) = 0;
  virtual void allocate_static_data_exc_vxc( int32_t nbf, int32_t nshells, bool do_vxc ) = 0;
  virtual void allocate_static_data_den( int32_t nbf, int32_t nshells ) = 0;
  virtual void allocate_static_data_exc_grad( int32_t nbf, int32_t nshells, int32_t natoms ) = 0;
  virtual void allocate_static_data_exx( int32_t nbf, int32_t nshells, size_t nshell_pairs, size_t nprim_pair_total, int32_t max_l ) = 0;
  virtual void allocate_static_data_exx_ek_screening( size_t ntasks, int32_t nbf, int32_t nshells, int nshell_pairs, int32_t max_l ) = 0;

  // Send persistent data from host to device
  virtual void send_static_data_weights( const Molecule& mol, const MolMeta& meta ) = 0;
  virtual void send_static_data_density_basis( const double* P, int32_t ldp, const BasisSet<double>& basis ) = 0;
  virtual void send_static_data_shell_pairs( const BasisSet<double>&, const ShellPairCollection<double>& ) = 0;
  virtual void send_static_data_exx_ek_screening( const double* V_max, int32_t ldv, const BasisSetMap&, const ShellPairCollection<double>& ) = 0;

  /// Zero out the density integrands in device memory
  virtual void zero_den_integrands() = 0;

  /// Zero out the EXC / VXC integrands in device memory
  virtual void zero_exc_vxc_integrands() = 0;

  /// Zero out the EXC Gradient integrands in device memory
  virtual void zero_exc_grad_integrands() = 0;

  /// Zero out the EXX integrands in device memory
  virtual void zero_exx_integrands() = 0;

  /// Zero out intermediates for EXX EK screening
  virtual void zero_exx_ek_screening_intermediates() = 0;

  /** Generate task batch to execute on device
   *
   *  Generate a batch of XC tasks to execute on the device and 
   *  populate device memory for said batch.
   *
   *  TODO: this will depend on the integrand, we should refactor this
   *  to only allocate what is needed
   *
   *  @param[in] basis_map  Basis set map instance for passed basis object
   *                        (TODO, this should probably persist to avoid clashes)
   *  @param[in] task_begin Start iterator for XC task queue
   *  @param[in] task_end   End iterator for XC task queue
   *
   *  @returns iterator to last XC task queue which was not kept in the
   *           allocated batch (!= task_end)
   */
  virtual host_task_iterator generate_buffers( integrator_term_tracker terms,
    const BasisSetMap& basis_map, host_task_iterator task_begin,
    host_task_iterator task_end ) = 0;

  /** Retreive EXC/VXC integrands from device memory
   *
   *  @param[out] EXC  Integrated XC energy (host) for XC task
   *  @param[out] N_EL Integrated # electrons (host) for XC queue (accuracy metric)
   *  @param[out[ VXC  Integrated XC potential (host) for XC queue
   */
  virtual void retrieve_exc_vxc_integrands( double* EXC, double* N_EL,
    double* VXC, int32_t ldvxc ) = 0;

  /** Retreive EXC Gradient integrands from device memory
   *
   *  @param[out] EXC_GRAD  Integrated XC Gradient (host) for XC task
   *  @param[out] N_EL      Integrated # electrons (host) for XC queue 
   */
  virtual void retrieve_exc_grad_integrands( double* EXC_GRAD, double* N_EL ) = 0;

  /** Retreive Density integrands from device memory
   *
   *  @param[out] N_EL      Integrated # electrons (host) for XC queue 
   */
  virtual void retrieve_den_integrands( double* N_EL ) = 0;


  virtual void retrieve_exx_integrands( double* K, int32_t ldk ) = 0;

  virtual void retrieve_exx_ek_max_bfn_sum( double* MBS, int32_t nt) = 0;


  virtual void copy_weights_to_tasks( host_task_iterator task_begin, host_task_iterator task_end ) = 0;
  virtual void populate_submat_maps ( size_t, host_task_iterator begin, host_task_iterator end, const BasisSetMap& ) = 0;

  virtual double* vxc_device_data() = 0;
  virtual double* exc_device_data() = 0;
  virtual double* nel_device_data() = 0;
  virtual double* exx_k_device_data() = 0;
  virtual device_queue queue() = 0;


};

}
