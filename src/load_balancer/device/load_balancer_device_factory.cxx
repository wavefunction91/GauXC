#include "load_balancer_impl.hpp"
#include "load_balancer_device_factory.hpp"

#ifdef GAUXC_ENABLE_CUDA
#include "cuda/replicated_cuda_load_balancer.hpp"
#endif

namespace GauXC {

std::shared_ptr<LoadBalancer> LoadBalancerDeviceFactory::get_shared_instance(
  std::string kernel_name, const RuntimeEnvironment& rt,
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis,
  size_t pv
) {

  std::transform(kernel_name.begin(), kernel_name.end(), 
    kernel_name.begin(), ::toupper );


  if( kernel_name == "DEFAULT" ) kernel_name = "REPLICATED";

  std::unique_ptr<detail::LoadBalancerImpl> ptr = nullptr;
  #ifdef GAUXC_ENABLE_CUDA
  if( kernel_name == "REPLICATED" ) {
    ptr = std::make_unique<detail::DeviceReplicatedLoadBalancer>(
      rt, mol, mg, basis, pv
    );
  }
  #endif

  if( ! ptr ) GAUXC_GENERIC_EXCEPTION("Load Balancer Kernel Not Recognized: " + kernel_name);

  return std::make_shared<LoadBalancer>(std::move(ptr));

}

}

