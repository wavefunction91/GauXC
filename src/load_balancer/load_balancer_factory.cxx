#include "load_balancer_impl.hpp"
#include "host/replicated_load_balancer.hpp"
#include "cuda/replicated_load_balancer.hpp"

#ifdef GAUXC_ENABLE_MPI
  #define GAUXC_MPI_ARG    MPI_Comm comm,
  #define GAUXC_MPI_PARAM  comm,
#else
  #define GAUXC_MPI_ARG   
  #define GAUXC_MPI_PARAM 
#endif

namespace GauXC {

LoadBalancerFactory::LoadBalancerFactory( ExecutionSpace ex, std::string kernel_name ) :
  ex_(ex), kernel_name_(kernel_name) { }

std::shared_ptr<LoadBalancer> LoadBalancerFactory::get_shared_instance(
  GAUXC_MPI_ARG
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
) {

  std::unique_ptr<detail::LoadBalancerImpl> ptr;
  switch(ex_) {

    case ExecutionSpace::Host:
      ptr = std::make_unique<detail::HostReplicatedLoadBalancer>(
        GAUXC_MPI_PARAM mol, mg, basis
      );
      break;

    default:
      throw std::runtime_error("Unrecognized LB space");

  }

  return std::make_shared<LoadBalancer>(std::move(ptr));

}

LoadBalancer LoadBalancerFactory::get_instance(
  GAUXC_MPI_ARG
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
) {

  auto ptr = get_shared_instance(GAUXC_MPI_PARAM mol,mg,basis);
  return LoadBalancer(std::move(*ptr));

}


}

