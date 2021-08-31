#pragma once
#include "replicated_host_load_balancer.hpp"

namespace GauXC  {
namespace detail {

struct PetiteHostReplicatedLoadBalancer : public HostReplicatedLoadBalancer {

  template <typename... Args>
  PetiteHostReplicatedLoadBalancer( Args&&... args ):
    HostReplicatedLoadBalancer( std::forward<Args>(args)... ) { }

  PetiteHostReplicatedLoadBalancer( const PetiteHostReplicatedLoadBalancer& );
  PetiteHostReplicatedLoadBalancer( PetiteHostReplicatedLoadBalancer&& ) noexcept;

  ~PetiteHostReplicatedLoadBalancer() noexcept;

  std::unique_ptr<LoadBalancerImpl> clone() const override final;

  std::pair< std::vector<int32_t>, size_t > micro_batch_screen(
    const BasisSet<double>&, const std::array<double,3>&,
    const std::array<double,3>& ) const override final;

};

}
}
