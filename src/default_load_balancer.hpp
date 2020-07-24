#pragma once

#include "load_balancer_impl.hpp"

namespace GauXC  {
namespace detail {

class DefaultLoadBalancer : public LoadBalancerImpl {

protected:

  using basis_type = BasisSet<double>;
  std::vector< XCTask > create_local_tasks_() const override;

public:

  DefaultLoadBalancer() = delete;
  template <typename... Args>
  DefaultLoadBalancer( Args&&... args ):
    LoadBalancerImpl( std::forward<Args>(args)... ) { }

  DefaultLoadBalancer( const DefaultLoadBalancer& );
  DefaultLoadBalancer( DefaultLoadBalancer&& ) noexcept;

  virtual ~DefaultLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override;

};

}
}
