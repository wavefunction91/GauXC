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
#include <gauxc/external/hdf5.h>
#include <gauxc/external/hdf5.hpp>

#include "c_status.hpp"
#include "c_molecule.hpp"
#include "c_basisset.hpp"
#include "c_matrix.hpp"
#include "hdf5_util.hpp"

namespace GauXC::C {
extern "C" {

/**
 * @brief Read a Molecule record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param mol Handle to the Molecule to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
void gauxc_molecule_read_hdf5_record(
  GauXCStatus* status,
  GauXCMolecule mol,
  const char* fname,
  const char* dset
) {
  detail::gauxc_status_init(status);
  try {
    read_hdf5_record( *detail::get_molecule_ptr(mol), std::string(fname), std::string(dset) );
  } catch(std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

/**
 * @brief Read a BasisSet record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param basis Handle to the BasisSet to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
void gauxc_basisset_read_hdf5_record(
  GauXCStatus* status,
  GauXCBasisSet basis,
  const char* fname,
  const char* dset
) {
  detail::gauxc_status_init(status);
  try {
    read_hdf5_record( *detail::get_basisset_ptr(basis), std::string(fname), std::string(dset) );
  } catch(std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

/**
 * @brief Read a CMatrix record from an HDF5 file.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the CMatrix to read into.
 * @param fname Name of the HDF5 file.
 * @param dset Name of the dataset within the HDF5 file.
 */
void gauxc_matrix_read_hdf5_record(
  GauXCStatus* status,
  GauXCMatrix matrix,
  const char* fname,
  const char* dset
) {
  detail::gauxc_status_init(status);
  detail::CMatrix& mat = *detail::get_matrix_ptr(matrix);
  try {
    HighFive::File file( std::string(fname), HighFive::File::ReadOnly );

    auto dataset = file.getDataSet(std::string(dset));
    auto dims = dataset.getDimensions();

    mat.resize( dims[0], dims[1] );

    dataset.read(mat.data());
  } catch(std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

} // extern "C"
} // namespace GauXC::C