#pragma once

#include "xc_device_data.hpp"

namespace GauXC {

struct XCDeviceStackData : public XCDeviceData {

  using XCDeviceData::host_task_type;
  using XCDeviceData::host_task_container;
  using XCDeviceData::host_task_iterator;

  size_t nshells = 0; ///< Number of shells allocated for static data
  size_t nbf     = 0; ///< Number of bfns allocated for static data
  size_t natoms  = 0; ///< Number of atoms allocated for static data
  
  void* device_ptr = nullptr; ///< Device buffer for all device allocations
  void* dynmem_ptr = nullptr; ///< Device buffer for dynamic allocations (mod static)
  size_t devmem_sz = 0;       ///< Length of device_ptr in bytes
  size_t dynmem_sz = 0;       ///< Length of dynmem_ptr in bytes 

  // Stack static data (not dynamically allocated for each task batch)

  Shell<double>* shells_device = nullptr; ///< Array of static basis shells (nshells)

  double* dmat_device   = nullptr; ///< Static density matrix storage (nbf,nbf)
  double* rab_device    = nullptr; ///< Static RAB matrix storage (*,natoms)
  double* coords_device = nullptr; ///< Static atomic positions (3 * natoms)

  double* exc_device     = nullptr; ///< EXC storage (1)
  double* nel_device     = nullptr; ///< N_EL storage (1)
  double* vxc_device     = nullptr; ///< VXC storage (nbf,nbf)
  double* acc_scr_device = nullptr; ///< Accumulaion scratch (1)


  // Stack dynamic data

  size_t total_npts_task_batch = 0; ///< Number of grid points in task batch

  double* points_device  = nullptr; ///< Grid points for task batch
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


  virtual ~XCDeviceStackData() noexcept = default;

  // Final overrides
  host_task_iterator generate_buffers( const BasisSetMap&,
    host_task_iterator, host_task_iterator) override final;
  void allocate_static_data( int32_t, int32_t, int32_t ) override final;
  void send_static_data( const double* P, int32_t ldp,
    const BasisSet<double>& basis, const Molecule& mol,
    const MolMeta& meta ) override final;
  void zero_integrands() override final;
  void retrieve_xc_integrands( double* EXC, double* N_EL,
    double* VXC, int32_t ldvxc ) override final;


  // New overridable APIs
  using device_buffer_t = std::tuple<void*, size_t>;
  virtual device_buffer_t alloc_pack_and_send( host_task_iterator begin, 
    host_task_iterator end, device_buffer_t buf, const BasisSetMap& );

  virtual size_t get_mem_req( const host_task_type&, const BasisSetMap& );


  // Device specific APIs
  virtual void allocate_device_buffer(int64_t sz = -1) = 0;
  virtual void master_queue_synchronize() = 0;
  virtual void copy_async_( size_t, const void*, void*, std::string ) = 0;
  virtual void set_zero_( size_t, void*, std::string );

  // Implementation specific APIs
  virtual void allocate_rab() = 0;
  virtual void send_rab(const MolMeta&) = 0;
  virtual size_t get_submat_chunk_size(int32_t,int32_t) = 0;
  virtual size_t get_static_mem_requirement() = 0;


  template <typename T>
  void copy_async( size_t sz, const T* src, T* dest, std::string msg ) {
    copy_async_( sz * sizeof(T), src, dest, msg );
  }

  template <typename T>
  void set_zero(size_t sz, T* data, std::string msg) {
    set_zero_(sz*sizeof(T), data, msg);
  }
};

}
