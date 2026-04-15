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
#include <gauxc/c/status.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif
 
enum GauXC_Type {
  GauXC_Type_Molecule = 1,
  GauXC_Type_BasisSet = 2,
  GauXC_Type_MolGrid  = 3,
  GauXC_Type_RuntimeEnvironment = 4,
  GauXC_Type_LoadBalancer = 5,
  GauXC_Type_LoadBalancerFactory = 6,
  GauXC_Type_MolecularWeights = 7,
  GauXC_Type_MolecularWeightsFactory = 8,
  GauXC_Type_Functional = 9,
  GauXC_Type_Integrator = 10,
};

/**
 * @brief GauXC C API Header for all GauXC objects.
 */
typedef struct GauXCHeader {
  enum GauXC_Type type; ///< Type of the GauXC object.
} GauXCHeader;

/**
 * @brief Delete a GauXC object.
 *
 * @param status Status object to capture any errors.
 * @param handle Address of the GauXC handle object pointer to delete
 * 
 * The function will delete the underlying GauXC object associated with the
 * handle and may set the handle to `NULL` on success.
 */
extern void gauxc_object_delete(
  GauXCStatus* status,
  void** handle
);

/**
 * @brief Delete all GauXC objects.
 *
 * @param status Status object to capture any errors.
 * @param handles Array of pointers to GauXC handle objects.
 * @param nhandles Number of handles in the array.
 * 
 * For each element, the underlying GauXC object will be deleted and the
 * corresponding handle pointer may be set to `NULL` on success.
 */
extern void gauxc_objects_delete(
  GauXCStatus* status,
  void** handles,
  size_t nhandles
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif