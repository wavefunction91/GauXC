/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/load_balancer.hpp>

namespace GauXC {

struct LoadBalancerDeviceFactory {

  static std::shared_ptr<LoadBalancer> get_shared_instance(
    std::string kernel_name, const RuntimeEnvironment& rt, 
    const Molecule& mol, const MolGrid& mg, const BasisSet<double>& basis
  );

};


}
