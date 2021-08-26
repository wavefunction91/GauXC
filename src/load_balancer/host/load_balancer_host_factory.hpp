#include "load_balancer_impl.hpp"

namespace GauXC {

struct LoadBalancerHostFactory {

  std::shared_ptr<LoadBalancer> get_shared_instance(
    #ifdef GAUXC_ENABLE_MPI
    MPI_Comm comm,
    #endif
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
  );

};


}
