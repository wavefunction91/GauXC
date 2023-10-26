/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_host_integrator.hpp>
#include "host/reference_replicated_xc_host_integrator.hpp"
#include "shell_batched_pgas_dist_xc_integrator.hpp"

namespace GauXC {
namespace detail {

template <typename ValueType>
class ShellBatchedPGASDistributedXCHostIntegrator : 
  public ShellBatchedPGASDistributedXCIntegrator<
    PGASDistributedXCHostIntegrator<ValueType>,
    ReferenceReplicatedXCHostIntegrator<ValueType>
  > {

  using base_type  = ShellBatchedPGASDistributedXCIntegrator<
    PGASDistributedXCHostIntegrator<ValueType>,
    ReferenceReplicatedXCHostIntegrator<ValueType>
  >;

public:

  template <typename... Args>
  ShellBatchedPGASDistributedXCHostIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ShellBatchedPGASDistributedXCHostIntegrator() noexcept;

};

extern template class ShellBatchedPGASDistributedXCHostIntegrator<double>;

}
}

