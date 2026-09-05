/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/c/xc_integrator.h>

#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>

namespace GauXC::detail {

static inline ReplicatedXCIntegratorImpl<double>* get_xc_integrator_ptr(C::GauXCIntegrator integrator) noexcept {
  return static_cast<ReplicatedXCIntegratorImpl<double>*>(integrator.ptr);
}

} // namespace GauXC::detail