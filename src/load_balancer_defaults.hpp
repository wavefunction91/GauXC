#include "default_load_balancer.hpp"

namespace GauXC {
namespace detail {

template <typename... Args>
std::unique_ptr<LoadBalancerImpl> make_default_load_balancer(Args&&... args) {
  return std::make_unique<DefaultLoadBalancer>( std::forward<Args>(args)... );
}

}
}
