#pragma once
#include "reduction_driver_impl.hpp"


namespace GauXC {

struct HostReductionDriver : public detail::ReductionDriverImpl {

  bool takes_host_memory() const override; 
  bool takes_device_memory() const override;

  virtual ~HostReductionDriver() noexcept;

  HostReductionDriver(const RuntimeEnvironment& rt);

};

}
