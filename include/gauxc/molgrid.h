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
#else
#include <stdint.h>
#endif
#include <gauxc/types.h>
#include <gauxc/status.h>
#include <gauxc/enums.h>
#include <gauxc/molecule.h>

#ifdef __cplusplus
extern "C" {
namespace GauXC::C {
#endif

/**
 * @brief GauXC C API MolGrid handle.
 */
typedef struct GauXCMolGrid {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the MolGrid instance.
} GauXCMolGrid;

/**
 * @brief Create a new MolGrid instance from atomic grids.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule.
 * @param pruning_scheme Pruning scheme to use.
 * @param batchsize Batch size for grid generation.
 * @param radial_quad Radial quadrature scheme to use.
 * @param grid_size Default atomic grid size to use.
 * @return Handle to the created MolGrid.
 */
extern GauXCMolGrid gauxc_molgrid_new_default(
    GauXCStatus* status,
    const GauXCMolecule mol,
    enum GauXC_PruningScheme pruning_scheme,
    int64_t batchsize,
    enum GauXC_RadialQuad radial_quad,
    enum GauXC_AtomicGridSizeDefault grid_size
);

/**
 * @brief Delete a MolGrid instance.
 * @param status Status object to capture any errors.
 * @param molgrid Handle to the MolGrid to delete.
 */
extern void gauxc_molgrid_delete(
    GauXCStatus* status,
    GauXCMolGrid* molgrid
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif