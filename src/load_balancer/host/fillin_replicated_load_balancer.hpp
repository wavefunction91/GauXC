#pragma once
#include "replicated_host_load_balancer.hpp"

namespace GauXC  {
namespace detail {

struct FillInHostReplicatedLoadBalancer : public HostReplicatedLoadBalancer {

  template <typename... Args>
  FillInHostReplicatedLoadBalancer( Args&&... args ):
    HostReplicatedLoadBalancer( std::forward<Args>(args)... ) { }

  FillInHostReplicatedLoadBalancer( const FillInHostReplicatedLoadBalancer& );
  FillInHostReplicatedLoadBalancer( FillInHostReplicatedLoadBalancer&& ) noexcept;

  ~FillInHostReplicatedLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override final;

  std::pair< std::vector<int32_t>, size_t > micro_batch_screen(
    const BasisSet<double>&, const std::array<double,3>&,
    const std::array<double,3>& ) const override final;

};

}
}
