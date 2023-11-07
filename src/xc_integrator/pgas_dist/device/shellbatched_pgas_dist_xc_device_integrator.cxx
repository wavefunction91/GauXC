/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "shellbatched_pgas_dist_xc_device_integrator.hpp"
#include "shell_batched_pgas_dist_xc_integrator_integrate_den.hpp"
#include "shell_batched_pgas_dist_xc_integrator_exc_vxc.hpp"
#include "shell_batched_pgas_dist_xc_integrator_exc_grad.hpp"
#include "shell_batched_pgas_dist_xc_integrator_exx.hpp"

namespace GauXC  {
namespace detail {

template <typename ValueType>
ShellBatchedPGASDistributedXCDeviceIntegrator<ValueType>::~ShellBatchedPGASDistributedXCDeviceIntegrator() noexcept = default;

template class ShellBatchedPGASDistributedXCDeviceIntegrator<double>;

}
}