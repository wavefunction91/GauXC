#include "load_balancer/host/replicated_load_balancer.hpp"
#include "load_balancer/cuda/replicated_load_balancer.hpp"

namespace GauXC {
namespace detail {

template <typename... Args>
std::unique_ptr<LoadBalancerImpl> make_default_load_balancer(Args&&... args) {
//#ifdef GAUXC_ENABLE_CUDA
//  return std::make_unique<DeviceReplicatedLoadBalancer>( std::forward<Args>(args)... );
//#else 
  return std::make_unique<HostReplicatedLoadBalancer>( std::forward<Args>(args)... );
//#endif
}

}
}
