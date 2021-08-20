#include <gauxc/oop_xc_integrator/replicated/replicated_xc_device_integrator.hpp>
#include "incore_replicated_xc_device_integrator.hpp"
#include "device/local_device_work_driver.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReplicatedXCDeviceIntegrator<ValueType>::~ReplicatedXCDeviceIntegrator() noexcept = default;

template class ReplicatedXCDeviceIntegrator<double>;


template <typename ValueType>
typename ReplicatedXCDeviceIntegratorFactory<ValueType>::ptr_return_t
  ReplicatedXCDeviceIntegratorFactory<ValueType>::make_integrator_impl(
    std::string integrator_kernel,
    std::shared_ptr<functional_type> func,
    std::shared_ptr<LoadBalancer> lb, 
    std::unique_ptr<LocalWorkDriver>&& lwd) {

  // Make sure that the LWD is a valid LocalDeviceWorkDriver
  if(not dynamic_cast<LocalDeviceWorkDriver*>(lwd.get())) {
    throw std::runtime_error("Passed LWD Not valid for Device ExSpace");
  }

  std::transform(integrator_kernel.begin(), integrator_kernel.end(), 
    integrator_kernel.begin(), ::toupper );

  if( integrator_kernel == "DEFAULT" ) integrator_kernel = "INCORE";

  if( integrator_kernel == "INCORE" )
    return std::make_unique<IncoreReplicatedXCDeviceIntegrator<ValueType>>(
      func, lb, std::move(lwd)
    );

  else
    throw std::runtime_error("INTEGRATOR KERNEL NOT RECOGNIZED");

  return nullptr;


}

template class ReplicatedXCDeviceIntegratorFactory<double>;


}
}

