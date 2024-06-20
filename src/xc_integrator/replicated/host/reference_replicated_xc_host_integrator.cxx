/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include "reference_replicated_xc_host_integrator_integrate_den.hpp"
#include "reference_replicated_xc_host_integrator_exc.hpp"
#include "reference_replicated_xc_host_integrator_exc_vxc.hpp"
#include "reference_replicated_xc_host_integrator_exc_vxc_neo.hpp"
#include "reference_replicated_xc_host_integrator_exc_grad.hpp"
#include "reference_replicated_xc_host_integrator_exx.hpp"
 
namespace GauXC::detail {

template <typename ValueType>
ReferenceReplicatedXCHostIntegrator<ValueType>::~ReferenceReplicatedXCHostIntegrator() noexcept = default;

template class ReferenceReplicatedXCHostIntegrator<double>;

}
