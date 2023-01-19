#include "device_reduction_driver.hpp"

namespace GauXC {

DeviceReductionDriver::DeviceReductionDriver(const RuntimeEnvironment& rt) :
  detail::ReductionDriverImpl(rt) { }


DeviceReductionDriver::~DeviceReductionDriver() noexcept = default;



bool DeviceReductionDriver::takes_host_memory() const {return false;}; 
bool DeviceReductionDriver::takes_device_memory() const {return true;};

}
