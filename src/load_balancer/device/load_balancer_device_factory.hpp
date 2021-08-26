#include <gauxc/load_balancer.hpp>

namespace GauXC {

struct LoadBalancerDeviceFactory {

  static std::shared_ptr<LoadBalancer> get_shared_instance(
    std::string kernel_name, 
    #ifdef GAUXC_ENABLE_MPI
    MPI_Comm comm,
    #endif
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
  );

};


}
