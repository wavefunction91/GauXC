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
#include <gauxc/status.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
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
  GauXC_Type_IntegratorFactory = 11,
  GauXC_Type_Matrix = 12
};

/**
 * @brief GauXC C API Header for all GauXC objects.
 */
typedef struct GauXCHeader {
  enum GauXC_Type type; ///< Type of the GauXC object.
} GauXCHeader;

/**
 * @brief Delete a GauXC object.
 * @param status Status object to capture any errors.
 * @param ptr Pointer to the GauXC object to delete.
 */
extern void gauxc_object_delete(
  GauXCStatus* status,
  void** ptr
);

/**
 * @brief Delete all GauXC objects.
 * @param status Status object to capture any errors.
 * @param ptrs Array of GauXCHeader objects.
 * @param nptrs Number of pointers in the array.
 */
extern void gauxc_objects_delete(
  GauXCStatus* status,
  void** ptrs,
  size_t nptrs
);

#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif