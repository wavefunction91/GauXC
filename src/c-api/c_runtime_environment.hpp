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

#include <gauxc/runtime_environment.h>
#include <gauxc/runtime_environment.hpp>

namespace GauXC::detail {
static inline RuntimeEnvironment* get_runtime_environment_ptr(C::GauXCRuntimeEnvironment env) noexcept {
#ifdef GAUXC_HAS_DEVICE
  void* ptr = env.device_ptr ? env.device_ptr : env.ptr;
#else
  void* ptr = env.ptr;
#endif
  return static_cast<RuntimeEnvironment*>(ptr);
}
#ifdef GAUXC_HAS_DEVICE
static inline DeviceRuntimeEnvironment* get_device_runtime_environment_ptr(C::GauXCRuntimeEnvironment env) noexcept {
  return static_cast<DeviceRuntimeEnvironment*>(env.device_ptr);
}
#endif
} // namespace GauXC::detail