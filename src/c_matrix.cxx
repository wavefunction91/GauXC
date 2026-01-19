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
#include <gauxc/util/c_matrix.hpp>
#include <gauxc/util/c_status.hpp>

namespace GauXC::C {
extern "C" {

GauXCMatrix gauxc_matrix_empty(
  GauXCStatus* status
) {
  status->code = 0;
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
  status->code = 0;
  GauXCMatrix matrix;
  matrix.hdr = GauXCHeader{GauXC_Type_Matrix};
  try { // can throw std::bad_alloc
    matrix.ptr = new detail::CMatrix( rows, cols );
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
    matrix.ptr = nullptr;
  }
  return matrix;
}

void gauxc_matrix_resize(
    GauXCStatus* status,
    const GauXCMatrix matrix,
    const size_t rows,
    const size_t cols
) {
  status->code = 0;
  try { // can throw std::bad_alloc
    detail::get_matrix_ptr(matrix)->resize( rows, cols );
  } catch (std::exception& e) {
    status->code = 1;
    status->message = detail::strdup(e.what());
  }
}

void gauxc_matrix_set_zero(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  status->code = 0;
  detail::get_matrix_ptr(matrix)->setZero();
}

size_t gauxc_matrix_rows(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  status->code = 0;
  return detail::get_matrix_ptr(matrix)->rows();
}

size_t gauxc_matrix_cols(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  status->code = 0;
  return detail::get_matrix_ptr(matrix)->cols();
}

detail::CMatrix::value_type* gauxc_matrix_data(
    GauXCStatus* status,
    const GauXCMatrix matrix
) {
  status->code = 0;
  return detail::get_matrix_ptr(matrix)->data();
}

void gauxc_matrix_delete(
  GauXCStatus* status,
  GauXCMatrix* matrix
) {
  status->code = 0;
  if(matrix == nullptr) return;
  if(matrix->ptr != nullptr)
    delete detail::get_matrix_ptr(*matrix);
  matrix->ptr = nullptr;
}

} // extern "C"
} // namespace GauXC::C