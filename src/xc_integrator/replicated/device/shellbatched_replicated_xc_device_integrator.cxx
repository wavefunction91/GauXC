/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "shellbatched_replicated_xc_device_integrator.hpp"
#include "shell_batched_replicated_xc_integrator_integrate_den.hpp"
#include "shell_batched_replicated_xc_integrator_exc_vxc.hpp"
#include "shell_batched_replicated_xc_integrator_exc_grad.hpp"
#include "shell_batched_replicated_xc_integrator_exx.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ShellBatchedReplicatedXCDeviceIntegrator<ValueType>::~ShellBatchedReplicatedXCDeviceIntegrator() noexcept = default;

template class ShellBatchedReplicatedXCDeviceIntegrator<double>;

}
}
