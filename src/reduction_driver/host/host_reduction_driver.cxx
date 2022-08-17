#include "host_reduction_driver.hpp"

namespace GauXC {

HostReductionDriver::HostReductionDriver(const RuntimeEnvironment& rt) :
  detail::ReductionDriverImpl(rt) { }


HostReductionDriver::~HostReductionDriver() noexcept = default;



bool HostReductionDriver::takes_host_memory() const {return true;}; 
bool HostReductionDriver::takes_device_memory() const {return false;};

}
