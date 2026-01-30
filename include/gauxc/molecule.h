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
#include <gauxc/types.h>
#include <gauxc/status.h>
#include <gauxc/atom.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
#endif

/**
 * @brief GauXC C API Molecule handle.
 */
typedef struct GauXCMolecule {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the Molecule instance.
} GauXCMolecule;

/**
 * @brief Create a new empty Molecule instance.
 * @param status Status object to capture any errors.
 * @return Handle to the created Molecule.
 */
extern GauXCMolecule gauxc_molecule_new(GauXCStatus* status);

/**
 * @brief Create a new Molecule instance from an array of Atoms.
 * @param status Status object to capture any errors.
 * @param atoms Pointer to an array of GauXCAtom.
 * @param natoms Number of atoms in the array.
 * @return Handle to the created Molecule.
 */
extern GauXCMolecule gauxc_molecule_new_from_atoms(GauXCStatus* status, GauXCAtom* atoms, size_t natoms );

/**
 * @brief Delete a Molecule instance.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule to delete.
 */
extern void gauxc_molecule_delete(GauXCStatus* status, GauXCMolecule* mol );

/**
 * @brief Get the number of atoms in the Molecule.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule.
 * @return Number of atoms in the Molecule.
 */
extern size_t gauxc_molecule_natoms(GauXCStatus* status, const GauXCMolecule mol );

/**
 * @brief Check if two Molecule instances are equal.
 * @param status Status object to capture any errors.
 * @param mol1 Handle to the first Molecule.
 * @param mol2 Handle to the second Molecule.
 * @return true if the Molecules are equal, false otherwise.
 */
extern bool gauxc_molecule_equal(
  GauXCStatus* status,
  const GauXCMolecule mol1,
  const GauXCMolecule mol2
);

#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif