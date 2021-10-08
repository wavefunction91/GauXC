#include "load_balancer_impl.hpp"
#include "load_balancer_host_factory.hpp"
#include "petite_replicated_load_balancer.hpp"
#include "fillin_replicated_load_balancer.hpp"

namespace GauXC {

std::shared_ptr<LoadBalancer> LoadBalancerHostFactory::get_shared_instance(
  std::string kernel_name, GAUXC_MPI_CODE(MPI_Comm comm,)
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis,
  size_t pv
) {

  std::transform(kernel_name.begin(), kernel_name.end(), 
    kernel_name.begin(), ::toupper );


  if( kernel_name == "DEFAULT" or kernel_name == "REPLICATED" ) 
    kernel_name = "REPLICATED-PETITE";

  std::unique_ptr<detail::LoadBalancerImpl> ptr = nullptr;
  if( kernel_name == "REPLICATED-PETITE" )
    ptr = std::make_unique<detail::PetiteHostReplicatedLoadBalancer>(
      GAUXC_MPI_CODE(comm,) mol, mg, basis, pv
    );

  if( kernel_name == "REPLICATED-FILLIN" )
    ptr = std::make_unique<detail::FillInHostReplicatedLoadBalancer>(
      GAUXC_MPI_CODE(comm,) mol, mg, basis, pv
    );

  if( ! ptr ) GAUXC_GENERIC_EXCEPTION("Load Balancer Kernel Not Recognized: " + kernel_name);

  return std::make_shared<LoadBalancer>(std::move(ptr));

}

}

