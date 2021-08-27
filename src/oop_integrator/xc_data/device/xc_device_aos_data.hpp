#pragma once

#include "xc_device_stack_data.hpp"
#include "xc_device_task.hpp"

namespace GauXC {

/// Base type for XCDeviceData instances that address task batches as AoS
struct XCDeviceAoSData : public XCDeviceStackData {

  size_t total_nbe_sq_task_batch   = 0; ///< Sum of nbe*nbe for task batch
  size_t total_nbe_npts_task_batch = 0; ///< Sum of npts*nbe for task batch
  size_t total_nshells_task_batch  = 0; ///< Sum of nshells for task batch
  size_t total_ncut_task_batch     = 0; ///< Sum of ncut for task batch
  size_t total_nblock_task_batch   = 0; ///< Sum of nblock for task batch

  // Collocation buffers
  struct aos_stack_data {
    double* bf_eval_device    = nullptr; 
      ///< Contiguous batch local collocation for task batch
    double* dbf_x_eval_device = nullptr; ///< Derivative of `bf_eval_device` wrt x
    double* dbf_y_eval_device = nullptr; ///< Derivative of `bf_eval_device` wrt y
    double* dbf_z_eval_device = nullptr; ///< Derivative of `bf_eval_device` wrt z

    // VXC Z Matrix
    double* zmat_vxc_lda_gga_device = nullptr;
      ///< Contiguous batch local Z matrix for LDA/GGA VXC for task batch

    // Scratch buffer
    double* nbe_scr_device = nullptr; ///< nbe*nbe scratch allocated for task batch

    // AoS Buffers
    size_t*  shell_list_device   = nullptr;
      ///< Contiguous batch local shell left for task batch
    size_t*  shell_offs_device   = nullptr;
      ///< Contiguous batch local shell offsets for task batch
    int32_t* submat_cut_device   = nullptr;
      ///< Contiguous batch local submatrix cuts for task batch
    int32_t* submat_block_device = nullptr;
      ///< Contiguous batch local submatrix blocking factors for task batch

    // Indirection
    XCDeviceTask* device_tasks = nullptr; ///< Task indirection in device memory

    inline void reset() { std::memset(this,0,sizeof(aos_stack_data)); }
  };

  std::vector<XCDeviceTask> host_device_tasks; ///< Task indirection in host memory
  aos_stack_data aos_stack;

  XCDeviceAoSData() = delete;
  inline XCDeviceAoSData( std::unique_ptr<DeviceBackend>&& ptr ) :
    XCDeviceStackData( std::move(ptr) ) { }

  // Make it polymorphic
  virtual ~XCDeviceAoSData() noexcept = default;

  // AoS Specific API
 
  /** Get L2 compatiable submatrix block size for a specified matrix dimension
   *
   *  @param[in] LDA Leading dimension large matrix which is extracted from
   *  @param[in] dev_id ID of device to query memory information
   *
   *  @returns Submatrix blocking factor
   */
  virtual size_t get_submat_chunk_size( int32_t LDA, int32_t dev_id ) = 0;

  // Overridable API overrides
  virtual size_t get_mem_req( integrator_term_tracker, const host_task_type&, const BasisSetMap&) override;
  //virtual device_buffer_t alloc_pack_and_send( integrator_term_tracker, host_task_iterator begin, 
  //  host_task_iterator end, device_buffer_t buf, const BasisSetMap&) override;
  virtual device_buffer_t allocate_dynamic_stack( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, device_buffer_t buf,
    const BasisSetMap& basis_map ) override;
  virtual void pack_and_send( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, 
    const BasisSetMap& basis_map ) override;

  inline virtual void 
    add_extra_to_indirection(integrator_term_tracker, 
                             std::vector<XCDeviceTask>&) { };

  virtual void reset_allocations() override;


};

}
