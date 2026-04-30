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

#include <gauxc/c/functional.h>

#include <gauxc/xc_integrator.hpp>

namespace GauXC::detail {
static inline functional_type* get_functional_ptr(C::GauXCFunctional mol) noexcept {
  return static_cast<functional_type*>(mol.ptr);
}
}