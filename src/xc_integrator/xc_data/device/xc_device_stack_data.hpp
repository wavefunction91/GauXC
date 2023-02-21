/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#pragma once

#include "xc_device_data.hpp"
#include "xc_device_shell_pair_soa.hpp"
#include "device/device_backend.hpp"
#include <cstring>
#include <gauxc/runtime_environment/fwd.hpp>

namespace GauXC {

// Collection of dimensions used in the XC integration
struct allocated_dims {
  size_t nshells      = 0; ///< Number of shells allocated for static data
  size_t nshell_pairs = 0; ///< Number of shell pairs allocated for static data
  size_t nbf          = 0; ///< Number of bfns allocated for static data
  size_t natoms       = 0; ///< Number of atoms allocated for static data
  size_t max_l        = 0; ///< Highest angular momentum value used
  size_t ntask_ek     = 0; ///< Number of total tasks allocated for static data (EK)
};

/// Base type for XCDeviceData instances that use stack data allocation.
struct XCDeviceStackData : public XCDeviceData {

  using XCDeviceData::host_task_type;
  using XCDeviceData::host_task_container;
  using XCDeviceData::host_task_iterator;

  allocated_dims global_dims; ///< Global dimensions for allocated data structures
  integrator_term_tracker allocated_terms;
  
  void* device_ptr = nullptr; ///< Device buffer for all device allocations
  void* dynmem_ptr = nullptr; ///< Device buffer for dynamic allocations (mod static)
  size_t devmem_sz = 0;       ///< Length of device_ptr in bytes
  size_t dynmem_sz = 0;       ///< Length of dynmem_ptr in bytes 

  // Stack static data (not dynamically allocated for each task batch)

  struct static_data {
    Shell<double>* shells_device = nullptr; ///< Array of static basis shells (nshells)
    ShellPair<double>* shell_pairs_device = nullptr;

    double* dmat_device   = nullptr; ///< Static density matrix storage (nbf,nbf)
    double* rab_device    = nullptr; ///< Static RAB matrix storage (*,natoms)
    double* coords_device = nullptr; ///< Static atomic positions (3 * natoms)

    double* exc_device     = nullptr;  ///< EXC storage (1)
    double* nel_device     = nullptr;  ///< N_EL storage (1)
    double* vxc_device     = nullptr;  ///< VXC storage (nbf,nbf)
    double* exx_k_device   = nullptr;  ///< EXX K storage (nbf,nbf)
    double* acc_scr_device = nullptr;  ///< Accumulaion scratch (1)
    double* exc_grad_device = nullptr; ///< EXC Gradient storage (3*natoms)

    double* vshell_max_device = nullptr;
    double* ek_bfn_max_device = nullptr;
    double* max_f_device      = nullptr;

    inline void reset() { std::memset( this, 0, sizeof(static_data) ); }
  };

  XCDeviceShellPairSoA shell_pair_soa;
  static_data static_stack;


  // Stack dynamic data

  size_t total_npts_task_batch = 0; ///< Number of grid points in task batch
  struct base_stack_data {

    //double* points_device  = nullptr; ///< Grid points for task batch
    double* points_x_device = nullptr;
    double* points_y_device = nullptr;
    double* points_z_device = nullptr;
    double* weights_device = nullptr; ///< Grid weights for task batch

    // U variables
    double* den_eval_device   = nullptr; ///< density for task batch
    double* den_x_eval_device = nullptr; ///< d/dx density for task batch
    double* den_y_eval_device = nullptr; ///< d/dy density for task batch
    double* den_z_eval_device = nullptr; ///< d/dz density for task batch

    // V variables / XC output
    double* gamma_eval_device  = nullptr; ///< gamma for task batch
    double* eps_eval_device    = nullptr; ///< XC energy density for task batch
    double* vrho_eval_device   = nullptr; ///< Rho XC derivative for task batch
    double* vgamma_eval_device = nullptr; ///< Gamma XC derivative for task batch

    inline void reset() { std::memset( this, 0, sizeof(base_stack_data) ); }
  };

  base_stack_data base_stack;

  /// Device backend instance to handle device specific execution
  const DeviceRuntimeEnvironment& runtime_;
  DeviceBackend* device_backend_ = nullptr;

  XCDeviceStackData() = delete; // No default ctor, must have device backend
  XCDeviceStackData( const DeviceRuntimeEnvironment& rt );

  virtual ~XCDeviceStackData() noexcept;

  // Final overrides
  host_task_iterator generate_buffers( integrator_term_tracker, const BasisSetMap&,
    host_task_iterator, host_task_iterator) override final;
  void allocate_static_data_weights( int32_t natoms ) override final;
  void allocate_static_data_exc_vxc( int32_t nbf, int32_t nshells ) override final;
  void allocate_static_data_den( int32_t nbf, int32_t nshells ) override final;
  void allocate_static_data_exc_grad( int32_t nbf, int32_t nshells, int32_t natoms ) override final;
  void allocate_static_data_exx( int32_t nbf, int32_t nshells, size_t nshell_pairs, int32_t max_l ) override final;
  void allocate_static_data_exx_ek_screening( size_t ntasks, int32_t nbf, int32_t nshells, int32_t max_l ) override final;
  void send_static_data_weights( const Molecule& mol, const MolMeta& meta ) override final;
  void send_static_data_density_basis( const double* P, int32_t ldp, 
    const BasisSet<double>& basis ) override final;
  void send_static_data_shell_pairs( const BasisSet<double>&, const ShellPairCollection<double>& ) 
    override final;
  void send_static_data_exx_vshell_max( const double* V_max, int32_t ldv ) override final;
  void zero_den_integrands() override final;
  void zero_exc_vxc_integrands() override final;
  void zero_exc_grad_integrands() override final;
  void zero_exx_integrands() override final;
  void zero_exx_ek_screening_intermediates() override final;
  void retrieve_exc_vxc_integrands( double* EXC, double* N_EL,
    double* VXC, int32_t ldvxc ) override final;
  void retrieve_exc_grad_integrands( double* EXC_GRAD, double* N_EL ) override final;
  void retrieve_den_integrands( double* N_EL ) override final;
  void retrieve_exx_integrands( double* K, int32_t ldk ) override final;
  void retrieve_exx_ek_approx_fmax( double* FMAX, int32_t ldF ) override final;
  void copy_weights_to_tasks( host_task_iterator task_begin, host_task_iterator task_end ) override final;

  double* vxc_device_data() override;
  double* exc_device_data() override;
  double* nel_device_data() override;
  double* exx_k_device_data() override;
  device_queue queue() override;


  virtual void reset_allocations() override;

  // New overridable APIs
  using device_buffer_t = std::tuple<void*, size_t>;

  /** Allocate and populate device memory for a given task batch
   *
   *  Overridable in devrived classes - derived classes should call
   *  this function explicitly to ensure that the correct information
   *  is allocated on the stack
   *
   *  @param[in] begin      Start iterator for task batch
   *  @param[in] end        End iterator for task batch
   *  @param[in] buf        Current state of dynamic memory stack
   *  @param[in] basis_map  Basis map instance for pass basis set 
   *                        (TODO: should persist internally)
   *
   *  @returns The state of the dynamic memory stack after allocating
   *           base information.
   */


  virtual device_buffer_t allocate_dynamic_stack( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, device_buffer_t buf );

  virtual void pack_and_send( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, 
    const BasisSetMap& basis_map );



  /** Obtain the memory requirement for an XC task
   *
   *  Overridable in devrived classes - derived classes should call
   *  this function explicitly to ensure that the correct information
   *  is allocated on the stack
   *
   *  @param[in] task       Task to obtain the memory requirement
   *
   *  @returns Memory requirement (bytes) for `task` in device memory
   */
  virtual size_t get_mem_req( integrator_term_tracker terms,
    const host_task_type& task );


  // Implementation specific APIs
  virtual size_t get_ldatoms()   = 0; ///< Stride of RAB in device memory
  virtual size_t get_rab_align() = 0; ///< Alignment of RAB in device memory
  virtual int get_points_per_subtask() = 0; ///< Number of points per subtask for OS kernels
  virtual size_t get_static_mem_requirement() = 0;
    ///< Static memory requirment for task batch which is independent of batch size

};

}
