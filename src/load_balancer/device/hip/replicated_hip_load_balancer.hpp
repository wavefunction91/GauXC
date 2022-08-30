#pragma once

#include "load_balancer_impl.hpp"

namespace GauXC  {
namespace detail {

class DeviceReplicatedLoadBalancer : public LoadBalancerImpl {

protected:

  using basis_type = BasisSet<double>;
  std::vector< XCTask > create_local_tasks_() const override;

public:

  DeviceReplicatedLoadBalancer() = delete;
  template <typename... Args>
  DeviceReplicatedLoadBalancer( Args&&... args ):
    LoadBalancerImpl( std::forward<Args>(args)... ) { }

  DeviceReplicatedLoadBalancer( const DeviceReplicatedLoadBalancer& );
  DeviceReplicatedLoadBalancer( DeviceReplicatedLoadBalancer&& ) noexcept;

  virtual ~DeviceReplicatedLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override;

};

}
}
