#pragma once

#include "xc_device_stack_data.hpp"
#include "xc_device_task.hpp"

namespace GauXC {

struct XCDeviceAoSData : public XCDeviceStackData {

  size_t total_nbe_sq_task_batch   = 0;
  size_t total_nbe_npts_task_batch = 0;
  size_t total_nshells_task_batch  = 0;
  size_t total_ncut_task_batch     = 0;
  size_t total_nblock_task_batch   = 0;

  // Collocation buffers
  double* bf_eval_device    = nullptr;
  double* dbf_x_eval_device = nullptr;
  double* dbf_y_eval_device = nullptr;
  double* dbf_z_eval_device = nullptr;

  // VXC Z Matrix
  double* zmat_vxc_lda_gga_device = nullptr;

  // Scratch buffer
  double* nbe_scr_device = nullptr;

  // AoS Buffers
  size_t*  shell_list_device   = nullptr;
  size_t*  shell_offs_device   = nullptr;
  int32_t* submat_cut_device   = nullptr;
  int32_t* submat_block_device = nullptr;

  // Indirection
  XCDeviceTask* device_tasks = nullptr;
  std::vector<XCDeviceTask> host_device_tasks;

  // Make it polymorphic
  virtual ~XCDeviceAoSData() noexcept = default;

  // API overrides
  virtual size_t get_mem_req( const host_task_type&, const BasisSetMap&) override;
  virtual device_buffer_t alloc_pack_and_send( host_task_iterator begin, 
    host_task_iterator end, device_buffer_t buf, const BasisSetMap&) override;

};

}