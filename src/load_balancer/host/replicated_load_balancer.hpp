#pragma once

#include "load_balancer_impl.hpp"

namespace GauXC  {
namespace detail {

class HostReplicatedLoadBalancer : public LoadBalancerImpl {

protected:

  using basis_type = BasisSet<double>;
  std::vector< XCTask > create_local_tasks_() const override;

public:

  HostReplicatedLoadBalancer() = delete;
  template <typename... Args>
  HostReplicatedLoadBalancer( Args&&... args ):
    LoadBalancerImpl( std::forward<Args>(args)... ) { }

  HostReplicatedLoadBalancer( const HostReplicatedLoadBalancer& );
  HostReplicatedLoadBalancer( HostReplicatedLoadBalancer&& ) noexcept;

  virtual ~HostReplicatedLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override;

};

}
}
