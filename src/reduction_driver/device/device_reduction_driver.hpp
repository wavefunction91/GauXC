#pragma once
#include "reduction_driver_impl.hpp"


namespace GauXC {

struct DeviceReductionDriver : public detail::ReductionDriverImpl {

  bool takes_host_memory() const override; 
  bool takes_device_memory() const override;

  virtual ~DeviceReductionDriver() noexcept;

  DeviceReductionDriver(GAUXC_MPI_CODE(MPI_Comm comm));

};

}
