/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "reference_replicated_xc_host_integrator_integrate_den.hpp"
#include "reference_replicated_xc_host_integrator_exc_vxc.hpp"
#include "reference_replicated_xc_host_integrator_exc_grad.hpp"
#include "reference_replicated_xc_host_integrator_exx.hpp"
 
namespace GauXC  {
namespace detail {

template <typename ValueType>
ReferenceReplicatedXCHostIntegrator<ValueType>::~ReferenceReplicatedXCHostIntegrator() noexcept = default;

template class ReferenceReplicatedXCHostIntegrator<double>;

}
}
