#include "host_reduction_driver.hpp"

namespace GauXC {

HostReductionDriver::HostReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm)) :
  detail::ReductionDriverImpl(GAUXC_MPI_CODE(comm)) { }


HostReductionDriver::~HostReductionDriver() noexcept = default;



bool HostReductionDriver::takes_host_memory() const {return true;}; 
bool HostReductionDriver::takes_device_memory() const {return false;};

}
