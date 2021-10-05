#include "load_balancer_impl.hpp"
#include "host/load_balancer_host_factory.hpp"

#ifdef GAUXC_ENABLE_DEVICE
#include "device/load_balancer_device_factory.hpp"
#endif

namespace GauXC {

LoadBalancerFactory::LoadBalancerFactory( ExecutionSpace ex, std::string kernel_name ) :
  ex_(ex), kernel_name_(kernel_name) { }

std::shared_ptr<LoadBalancer> LoadBalancerFactory::get_shared_instance(
  GAUXC_MPI_CODE(MPI_Comm comm,)
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis,
  size_t pad_value
) {

  switch(ex_) {
    case ExecutionSpace::Host:
      using host_factory = LoadBalancerHostFactory;
      return host_factory::get_shared_instance(kernel_name_,
        GAUXC_MPI_CODE(comm,) mol, mg, basis, pad_value );
    #ifdef GAUXC_ENABLE_DEVICE
    case ExecutionSpace::Device:
      using device_factory = LoadBalancerDeviceFactory;
      return device_factory::get_shared_instance(kernel_name_,
        GAUXC_MPI_CODE(comm,) mol, mg, basis, pad_value );
    #endif
    default:
      throw std::runtime_error("Unrecognized LB space");
   }


}

LoadBalancer LoadBalancerFactory::get_instance(
  GAUXC_MPI_CODE(MPI_Comm comm,)
  const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis,
  size_t pad_value
) {

  auto ptr = get_shared_instance(GAUXC_MPI_CODE(comm,) mol,mg,basis, pad_value);
  return LoadBalancer(std::move(*ptr));

}


}

