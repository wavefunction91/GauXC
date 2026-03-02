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
#include <gauxc/c/runtime_environment.h>

#include <gauxc/runtime_environment.hpp>

#include "c_runtime_environment.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

GauXCRuntimeEnvironment gauxc_runtime_environment_new( 
  GauXCStatus* status
  GAUXC_MPI_CODE(, MPI_Comm comm) 
) {
  detail::gauxc_status_init(status);
  GauXCRuntimeEnvironment env{};
  env.ptr = nullptr;
#ifdef GAUXC_HAS_DEVICE
  env.device_ptr = nullptr;
#endif

  try {
    env.ptr = new RuntimeEnvironment(GAUXC_MPI_CODE(comm));
    env.hdr = GauXCHeader{GauXC_Type_RuntimeEnvironment};
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return env;
}

void gauxc_runtime_environment_delete(GauXCStatus* status, GauXCRuntimeEnvironment* env) {
  detail::gauxc_status_init(status);
  if (env == nullptr) return;
  if (env->ptr != nullptr)
    delete detail::get_runtime_environment_ptr(*env);
  env->ptr = nullptr;
#ifdef GAUXC_HAS_DEVICE
  if (env->device_ptr != nullptr)
    delete detail::get_device_runtime_environment_ptr(*env);
  env->device_ptr = nullptr;
#endif
}

int gauxc_runtime_environment_comm_rank(GauXCStatus* status, const GauXCRuntimeEnvironment env) {
  detail::gauxc_status_init(status);
  return detail::get_runtime_environment_ptr(env)->comm_rank();
}

int gauxc_runtime_environment_comm_size(GauXCStatus* status, const GauXCRuntimeEnvironment env) {
  detail::gauxc_status_init(status);
  return detail::get_runtime_environment_ptr(env)->comm_size();
}

#ifdef GAUXC_HAS_DEVICE
GauXCRuntimeEnvironment gauxc_device_runtime_environment_new( 
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  double fill_fraction
) {
  detail::gauxc_status_init(status);
  GauXCRuntimeEnvironment env{};
  env.hdr = GauXCHeader{GauXC_Type_RuntimeEnvironment};
  env.ptr = nullptr;
  env.device_ptr = nullptr;

  try {
    env.device_ptr = new DeviceRuntimeEnvironment(
      GAUXC_MPI_CODE(comm,) fill_fraction
    );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return env;
}

GauXCRuntimeEnvironment gauxc_device_runtime_environment_new_mem( 
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  void* mem,
  size_t mem_sz
) {
  detail::gauxc_status_init(status);
  GauXCRuntimeEnvironment env{};
  env.hdr = GauXCHeader{GauXC_Type_RuntimeEnvironment};
  env.ptr = nullptr;
  env.device_ptr = nullptr;

  try {
    env.device_ptr = new DeviceRuntimeEnvironment(
      GAUXC_MPI_CODE(comm,) mem, mem_sz
    );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return env;
}
#endif

} // extern "C"
} // namespace GauXC::C