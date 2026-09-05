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

#ifdef __cplusplus
#include <cstddef>
#else
#include <stddef.h>
#endif
#include <gauxc/c/types.h>
#include <gauxc/c/status.h>
#include <gauxc/c/mpi.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API RuntimeEnvironment handle.
 */
typedef struct GauXCRuntimeEnvironment {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the RuntimeEnvironment instance.
#ifdef GAUXC_HAS_DEVICE
  void* device_ptr; ///< Pointer to the DeviceRuntimeEnvironment instance (if applicable).
#endif
} GauXCRuntimeEnvironment;

/**
 * @brief Create a new RuntimeEnvironment instance.
 * @param status Status object to capture any errors.
 * @param comm MPI Communicator (if applicable).
 * @return Handle to the created RuntimeEnvironment.
 */
extern GauXCRuntimeEnvironment gauxc_runtime_environment_new(
  GauXCStatus* status
  GAUXC_MPI_CODE(, MPI_Comm comm) 
);

/**
 * @brief Delete a RuntimeEnvironment instance.
 * @param status Status object to capture any errors.
 * @param env Handle to the RuntimeEnvironment to delete.
 */
extern void gauxc_runtime_environment_delete( 
  GauXCStatus* status,
  GauXCRuntimeEnvironment* env 
);

/**
 * @brief Get the rank of the current process in the RuntimeEnvironment's communicator.
 * @param status Status object to capture any errors.
 * @param env Handle to the RuntimeEnvironment.
 * @return Rank of the current process.
 */
extern int gauxc_runtime_environment_comm_rank( 
  GauXCStatus* status,
  const GauXCRuntimeEnvironment env 
);
/**
 * @brief Get the size of the RuntimeEnvironment's communicator.
 * @param status Status object to capture any errors.
 * @param env Handle to the RuntimeEnvironment.
 * @return Size of the communicator.
 */
extern int gauxc_runtime_environment_comm_size(
  GauXCStatus* status,
  const GauXCRuntimeEnvironment env
);

#ifdef GAUXC_HAS_DEVICE
/**
 * @brief Create new DeviceRuntimeEnvironment instance.
 * @param status Status object to capture any errors.
 * @param comm MPI Communicator (if applicable).
 * @param fill_fraction Fraction of device memory to use.
 * @return Handle to the created DeviceRuntimeEnvironment.
 */
extern GauXCRuntimeEnvironment gauxc_device_runtime_environment_new(
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  double fill_fraction
);

/**
 * @brief Create new DeviceRuntimeEnvironment instance.
 * @param status Status object to capture any errors.
 * @param comm MPI Communicator (if applicable).
 * @param mem Pointer to preallocated device memory.
 * @param mem_sz Size of preallocated device memory.
 * @return Handle to the created DeviceRuntimeEnvironment.
 */
extern GauXCRuntimeEnvironment gauxc_device_runtime_environment_new_mem(
  GauXCStatus* status,
  GAUXC_MPI_CODE(MPI_Comm comm,)
  void* mem, 
  size_t mem_sz
);
#endif

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif