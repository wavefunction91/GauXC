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
#include <gauxc/gauxc_config.h>
#ifdef GAUXC_HAS_HDF5
#include <gauxc/status.h>
#include <gauxc/basisset.h>
#include <gauxc/molecule.h>
#include <gauxc/matrix.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
#endif

/**
 * @brief Write a Molecule record to an HDF5 file.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule to write.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_molecule_write_hdf5_record(
  GauXCStatus* status,
  GauXCMolecule mol,
  const char* fname,
  const char* dset
);

/**
 * @brief Write a BasisSet record to an HDF5 file.
 * @param status Status object to capture any errors.
 * @param basis Handle to the BasisSet to write.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_basisset_write_hdf5_record(
  GauXCStatus* status,
  GauXCBasisSet basis,
  const char* fname,
  const char* dset
);

/**
 * @brief Write a CMatrix record to an HDF5 file.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the CMatrix to write.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_matrix_write_hdf5_record(
  GauXCStatus* status,
  GauXCMatrix matrix,
  const char* fname,
  const char* dset
);

/**
 * @brief Read a Molecule record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_molecule_read_hdf5_record(
  GauXCStatus* status,
  GauXCMolecule mol,
  const char* fname,
  const char* dset
);

/**
 * @brief Read a BasisSet record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param basis Handle to the BasisSet to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_basisset_read_hdf5_record(
  GauXCStatus* status,
  GauXCBasisSet basis,
  const char* fname,
  const char* dset
);

/**
 * @brief Read a CMatrix record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the CMatrix to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
extern void gauxc_matrix_read_hdf5_record(
  GauXCStatus* status,
  GauXCMatrix matrix,
  const char* fname,
  const char* dset
);

#ifdef __cplusplus
} // extern "C"
} // namespace GauXC::C
#endif
#endif