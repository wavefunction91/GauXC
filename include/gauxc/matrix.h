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
#include <cstddef>
#else
#include <stdint.h>
#include <stddef.h>
#endif
#include <gauxc/types.h>
#include <gauxc/status.h>

#ifdef __cplusplus
namespace GauXC::C {
extern "C" {
#endif

/**
 * @brief GauXC C API Matrix handle.
 */
typedef struct GauXCMatrix {
  GauXCHeader hdr; ///< Header for internal use.
  void* ptr; ///< Pointer to the Matrix instance.
} GauXCMatrix;

/**
 * @brief Create a new Matrix instance.
 * @param status Status object to capture any errors.
 * @return Handle to the newly created Matrix.
 */
extern GauXCMatrix gauxc_matrix_empty(
    GauXCStatus* status
);

/**
 * @brief Create a new Matrix instance.
 * @param status Status object to capture any errors.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @return Handle to the newly created Matrix.
 */
extern GauXCMatrix gauxc_matrix_new(
    GauXCStatus* status,
    const size_t rows,
    const size_t cols
);

/**
 * @brief Resize an existing Matrix instance.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix to resize.
 * @param rows New number of rows in the matrix.
 * @param cols New number of columns in the matrix.
 */
extern void gauxc_matrix_resize(
    GauXCStatus* status,
    const GauXCMatrix matrix,
    const size_t rows,
    const size_t cols
);

/**
 * @brief Set all elements of the Matrix to zero.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix to set to zero.
 */
extern void gauxc_matrix_set_zero(
    GauXCStatus* status,
    const GauXCMatrix matrix
);

/**
 * @brief Get the number of rows in the Matrix.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix.
 * @return Number of rows in the Matrix.
 */
extern size_t gauxc_matrix_rows(
    GauXCStatus* status,
    const GauXCMatrix matrix
);

/**
 * @brief Get the number of columns in the Matrix.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix.
 * @return Number of columns in the Matrix.
 */
extern size_t gauxc_matrix_cols(
    GauXCStatus* status,
    const GauXCMatrix matrix
);

/**
 * @brief Get a pointer to the underlying data of the Matrix.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix.
 * @return Pointer to the underlying data of the Matrix.
 */
extern double* gauxc_matrix_data(
    GauXCStatus* status,
    const GauXCMatrix matrix
);

/**
 * @brief Delete a Matrix instance.
 * @param status Status object to capture any errors.
 * @param matrix Handle to the Matrix to delete.
 */
extern void gauxc_matrix_delete(
    GauXCStatus* status,
    GauXCMatrix* matrix
);

#ifdef __cplusplus
} // namespace GauXC::C
} // extern "C"
#endif