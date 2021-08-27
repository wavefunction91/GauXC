#include "device_reduction_driver.hpp"

namespace GauXC {

DeviceReductionDriver::DeviceReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm)) :
  detail::ReductionDriverImpl(GAUXC_MPI_CODE(comm)) { }


DeviceReductionDriver::~DeviceReductionDriver() noexcept = default;



bool DeviceReductionDriver::takes_host_memory() const {return false;}; 
bool DeviceReductionDriver::takes_device_memory() const {return true;};

}
