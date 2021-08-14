#include <gauxc/oop_xc_integrator/replicated/replicated_xc_host_integrator.hpp>
#include "reference_replicated_xc_host_integrator.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReplicatedXCHostIntegrator<ValueType>::~ReplicatedXCHostIntegrator() noexcept = default;

template class ReplicatedXCHostIntegrator<double>;


template <typename ValueType>
typename ReplicatedXCHostIntegratorFactory<ValueType>::ptr_return_t
  ReplicatedXCHostIntegratorFactory<ValueType>::make_integrator_impl(
    std::string integrator_kernel,
    std::shared_ptr<functional_type> func,
    std::shared_ptr<LoadBalancer> lb, 
    std::unique_ptr<LocalWorkDriver>&& lwd) {

  // Make sure that the LWD is a valid LocalHostWorkDriver
  if(not dynamic_cast<LocalHostWorkDriver*>(lwd.get())) {
    throw std::runtime_error("Passed LWD Not valid for Host ExSpace");
  }

  std::transform(integrator_kernel.begin(), integrator_kernel.end(), 
    integrator_kernel.begin(), ::toupper );

  if( integrator_kernel == "DEFAULT" ) integrator_kernel = "REFERENCE";

  if( integrator_kernel == "REFERENCE" )
    return std::make_unique<ReferenceReplicatedXCHostIntegrator<ValueType>>(
      func, lb, std::move(lwd)
    );

  else
    throw std::runtime_error("INTEGRATOR KERNEL NOT RECOGNIZED");

  return nullptr;


}

template class ReplicatedXCHostIntegratorFactory<double>;


}
}

