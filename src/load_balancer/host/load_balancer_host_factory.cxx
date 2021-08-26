#include "load_balancer_impl.hpp"
#include "load_balancer_host_factory.hpp"
#include "replicated_host_load_balancer.hpp"

#ifdef GAUXC_ENABLE_MPI
  #define GAUXC_MPI_ARG    MPI_Comm comm,
  #define GAUXC_MPI_PARAM  comm,
#else
  #define GAUXC_MPI_ARG   
  #define GAUXC_MPI_PARAM 
#endif

namespace GauXC {

std::shared_ptr<LoadBalancer> LoadBalancerHostFactory::get_shared_instance(
  std::string kernel_name, GAUXC_MPI_ARG
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
) {

  std::transform(kernel_name.begin(), kernel_name.end(), 
    kernel_name.begin(), ::toupper );


  if( kernel_name == "DEFAULT" ) kernel_name = "REPLICATED";

  std::unique_ptr<detail::LoadBalancerImpl> ptr = nullptr;
  if( kernel_name == "REPLICATED" )
    ptr = std::make_unique<detail::HostReplicatedLoadBalancer>(
      GAUXC_MPI_PARAM mol, mg, basis
    );

  if( ! ptr ) throw std::runtime_error("Load Balancer Kernel Not Recognized");

  return std::make_shared<LoadBalancer>(std::move(ptr));

}

}

