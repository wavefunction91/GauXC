#pragma once
#include "device/xc_device_aos_data.hpp"

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

  scheme1_data scheme1_stack;

  virtual ~Scheme1DataBase() noexcept;
  Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr);

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
