#pragma once
#include "device/xc_device_aos_data.hpp"

namespace GauXC {

struct Scheme1DataBase : public XCDeviceAoSData {

  using base_type = XCDeviceAoSData;
  using base_type::host_task_type;
  using base_type::device_buffer_t;

  double*  dist_scratch_device = nullptr;
  double*  dist_nearest_device = nullptr;
  int32_t* iparent_device      = nullptr;

  virtual ~Scheme1DataBase() noexcept;
  Scheme1DataBase(std::unique_ptr<DeviceBackend>&& ptr, bool batch_l3_blas = true);

  // Final overrides
  device_buffer_t add_extra_to_indirection(std::vector<XCDeviceTask>&, 
    device_buffer_t ) override final;

  // Overrideable API's
  virtual size_t get_mem_req( const host_task_type&, const BasisSetMap&) override;
  virtual size_t get_static_mem_requirement() override; 
};

}
