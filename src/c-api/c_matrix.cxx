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

#include <gauxc/matrix.h>

#include "c_matrix.hpp"
#include "c_status.hpp"

namespace GauXC::C {
extern "C" {

GauXCMatrix gauxc_matrix_empty(
  GauXCStatus* status
) {
  detail::gauxc_status_init(status);
  GauXCMatrix matrix;
  matrix.hdr = GauXCHeader{GauXC_Type_Matrix};
  matrix.ptr = new detail::CMatrix();

  return matrix;
}

GauXCMatrix gauxc_matrix_new(
  GauXCStatus* status,
  const size_t rows,
  const size_t cols
) {
  detail::gauxc_status_init(status);
  GauXCMatrix matrix;
  matrix.hdr = GauXCHeader{GauXC_Type_Matrix};
  try { // can throw std::bad_alloc
    matrix.ptr = new detail::CMatrix( rows, cols );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
  return matrix;
}

void gauxc_matrix_resize(
    GauXCStatus* status,
    const GauXCMatrix matrix,
    const size_t rows,
    const size_t cols
) {
  detail::gauxc_status_init(status);
  try { // can throw std::bad_alloc
    detail::get_matrix_ptr(matrix)->resize( rows, cols );
  } catch (std::exception& e) {
    detail::gauxc_status_handle(status, 1, e.what());
  }
}

void gauxc_matrix_set_zero(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  detail::gauxc_status_init(status);
  detail::get_matrix_ptr(matrix)->setZero();
}

size_t gauxc_matrix_rows(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  detail::gauxc_status_init(status);
  return detail::get_matrix_ptr(matrix)->rows();
}

size_t gauxc_matrix_cols(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  detail::gauxc_status_init(status);
  return detail::get_matrix_ptr(matrix)->cols();
}

detail::CMatrix::value_type* gauxc_matrix_data(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  detail::gauxc_status_init(status);
  return detail::get_matrix_ptr(matrix)->data();
}

void gauxc_matrix_delete(
  GauXCStatus* status,
  GauXCMatrix* matrix
) {
  detail::gauxc_status_init(status);
  if(matrix == nullptr) return;
  if(matrix->ptr != nullptr)
    delete detail::get_matrix_ptr(*matrix);
  matrix->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C