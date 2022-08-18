#pragma once
#include "device/xc_device_aos_data.hpp"
#include "device/common/shell_to_task.hpp"

namespace GauXC {

struct Scheme1DataBase : public XCDeviceAoSData {

  using base_type = XCDeviceAoSData;
  using base_type::host_task_type;
  using base_type::device_buffer_t;

  struct scheme1_data {
    double*  dist_scratch_device = nullptr;
    double*  dist_nearest_device = nullptr;
    int32_t* iparent_device      = nullptr;
    
    inline void reset(){ std::memset(this,0,sizeof(scheme1_data)); }
  };

  struct collocation_data {
    size_t*  shell_list_device   = nullptr;
      ///< Contiguous batch local shell left for task batch
    size_t*  shell_offs_device   = nullptr;
      ///< Contiguous batch local shell offsets for task batch
    inline void reset(){ std::memset(this,0,sizeof(collocation_data)); }
  };

  struct shell_to_task_data {
    ShellToTaskDevice* shell_to_task_device;

    int32_t* shell_to_task_idx_device = nullptr;
    int32_t* shell_to_task_off_device = nullptr;

    inline void reset(){ std::memset(this,0,sizeof(shell_to_task_data)); }
  };

  size_t total_nshells_task_batch  = 0; ///< Sum of nshells for task batch
  scheme1_data       scheme1_stack;
  collocation_data   collocation_stack;
  shell_to_task_data shell_to_task_stack;
  std::vector<AngularMomentumShellToTaskBatch> l_batched_shell_to_task;

  virtual ~Scheme1DataBase() noexcept;
  //Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr);
  Scheme1DataBase(const DeviceRuntimeEnvironment& rt);

  // Final overrides
  void add_extra_to_indirection(integrator_term_tracker, 
    std::vector<XCDeviceTask>& ) override final;

  // Overrideable API's
  virtual size_t get_mem_req( integrator_term_tracker, 
    const host_task_type&) override;
  virtual size_t get_static_mem_requirement() override; 
  virtual void reset_allocations() override;
  virtual device_buffer_t allocate_dynamic_stack( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, device_buffer_t buf )
    override;
  virtual void pack_and_send( integrator_term_tracker terms,
    host_task_iterator begin, host_task_iterator end, 
    const BasisSetMap& basis_map ) override;
};

}
