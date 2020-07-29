#pragma once

#include "load_balancer_impl.hpp"

namespace GauXC  {
namespace detail {

class ReplicatedLoadBalancer : public LoadBalancerImpl {

protected:

  using basis_type = BasisSet<double>;
  std::vector< XCTask > create_local_tasks_() const override;

public:

  ReplicatedLoadBalancer() = delete;
  template <typename... Args>
  ReplicatedLoadBalancer( Args&&... args ):
    LoadBalancerImpl( std::forward<Args>(args)... ) { }

  ReplicatedLoadBalancer( const ReplicatedLoadBalancer& );
  ReplicatedLoadBalancer( ReplicatedLoadBalancer&& ) noexcept;

  virtual ~ReplicatedLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override;

};

}
}
