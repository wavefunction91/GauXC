#pragma once

#include <gauxc/xc_integrator.hpp>

namespace GauXC::integrator::host {

void partition_weights_host(
  XCWeightAlg            weight_alg,
  const Molecule&        mol,
  const MolMeta&         meta,
  std::vector< XCTask >& tasks
);


}
