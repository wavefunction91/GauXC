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
#include <cstdint>
#include <cstdbool>
#include <cstddef>
#else
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#endif
#include <gauxc/c/types.h>
#include <gauxc/c/status.h>
#include <gauxc/c/shell.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API BasisSet handle.
 */
typedef struct GauXCBasisSet {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the BasisSet instance.
} GauXCBasisSet;

/**
 * @brief Create a new BasisSet instance.
 * @param status Status object to capture any errors.
 * @return Handle to the created BasisSet.
 */
extern GauXCBasisSet gauxc_basisset_new( GauXCStatus* status );

/**
 * @brief Create a new BasisSet instance from arrays of shells.
 * @param status Status object to capture any errors.
 * @param shells Pointer to an array of GauXCShell.
 * @param nshells Number of shells in the array.
 * @param normalize Whether to normalize the basis functions.
 * @return Handle to the created BasisSet.
 */
extern GauXCBasisSet gauxc_basisset_new_from_shells(
   GauXCStatus* status,
   GauXCShell* shells,
   size_t nshells,
   bool normalize
);

/**
 * @brief Delete a BasisSet instance.
 * @param status Status object to capture any errors.
 * @param basis Handle to the BasisSet to delete.
 */
extern void gauxc_basisset_delete( GauXCStatus* status, GauXCBasisSet* basis );

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif