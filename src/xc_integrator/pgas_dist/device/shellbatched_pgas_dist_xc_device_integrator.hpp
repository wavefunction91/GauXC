/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_device_integrator.hpp>
#include "device/incore_replicated_xc_device_integrator.hpp"
#include "shell_batched_pgas_dist_xc_integrator.hpp"

namespace GauXC {
namespace detail {

template <typename ValueType>
class ShellBatchedPGASDistributedXCDeviceIntegrator : 
  public ShellBatchedPGASDistributedXCIntegrator<
    PGASDistributedXCDeviceIntegrator<ValueType>,
    IncoreReplicatedXCDeviceIntegrator<ValueType>
  > {

  using base_type  = ShellBatchedPGASDistributedXCIntegrator<
    PGASDistributedXCDeviceIntegrator<ValueType>,
    IncoreReplicatedXCDeviceIntegrator<ValueType>
  >;

public:

  template <typename... Args>
  ShellBatchedPGASDistributedXCDeviceIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ShellBatchedPGASDistributedXCDeviceIntegrator() noexcept;

};

extern template class ShellBatchedPGASDistributedXCDeviceIntegrator<double>;

}
}
