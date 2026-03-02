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

#include <gauxc/c/molgrid.h>

#include <gauxc/molgrid.hpp>

namespace GauXC::detail {
static inline MolGrid* get_molgrid_ptr(C::GauXCMolGrid molgrid) noexcept {
  return static_cast<MolGrid*>(molgrid.ptr);
}
}