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

  virtual std::pair< std::vector<int32_t>, size_t > micro_batch_screen(
    const BasisSet<double>&, const std::array<double,3>&,
    const std::array<double,3>& ) const = 0;

};

}
}
