#pragma once
#include <gauxc/load_balancer.hpp>

namespace GauXC {

struct LoadBalancerHostFactory {

  static std::shared_ptr<LoadBalancer> get_shared_instance(
    std::string kernel_name, GAUXC_MPI_CODE(MPI_Comm comm,)
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis,
    size_t pv
  );

};


}
