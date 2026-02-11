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

#include <algorithm>
#include <cstddef>

namespace GauXC::detail {

class CMatrix {
public:
  using value_type = double;

  /**
   * @brief Default constructor creates a 0x0 matrix
   */
  CMatrix() noexcept : rows_(0), cols_(0), data_(nullptr) {}

  /**
   * @brief Construct a new CMatrix object
   * @param rows Number of rows
   * @param cols Number of columns
   */
  CMatrix(size_t rows, size_t cols)
      : rows_(rows), cols_(cols), data_(new value_type[rows * cols]()) {}

  /**
   * @brief Destroy the CMatrix object
   */
  ~CMatrix() noexcept {
    delete[] data_;  // delete[] on nullptr is safe
  }

  /**
   * @brief Copy constructor
   * @param other CMatrix to copy from
   */
  CMatrix(const CMatrix& other)
      : rows_(other.rows_), cols_(other.cols_), data_(nullptr) {
    if (other.data_ && rows_ * cols_ > 0) {
      data_ = new value_type[rows_ * cols_];
      std::copy(other.data_, other.data_ + rows_ * cols_, data_);
    }
  }

  /**
   * @brief Copy assignment operator
   * @param other CMatrix to copy from
   * @return CMatrix& Reference to this
   */
  CMatrix& operator=(const CMatrix& other) {
    if (this != &other) {
      // Always clean up if sizes differ OR if other has no data
      if (rows_ * cols_ != other.rows_ * other.cols_ || !other.data_) {
        delete[] data_;
        data_ = nullptr;
      }
      rows_ = other.rows_;
      cols_ = other.cols_;
      if (other.data_ && rows_ * cols_ > 0) {
        if (!data_) {
          data_ = new value_type[rows_ * cols_];
        }
        std::copy(other.data_, other.data_ + rows_ * cols_, data_);
      }
    }
    return *this;
  }

  /**
   * @brief Move constructor
   * @param other CMatrix to move from
   */
  CMatrix(CMatrix&& other) noexcept
      : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {
    other.rows_ = 0;
    other.cols_ = 0;
    other.data_ = nullptr;
  }

  // Move assignment
  /**
   * @brief Move assignment operator
   * @param other CMatrix to move from
   * @return CMatrix& Reference to this
   */
  CMatrix& operator=(CMatrix&& other) noexcept {
    if (this != &other) {
      delete[] data_;
      rows_ = other.rows_;
      cols_ = other.cols_;
      data_ = other.data_;
      other.rows_ = 0;
      other.cols_ = 0;
      other.data_ = nullptr;
    }
    return *this;
  }

  /**
   * @brief Get pointer to raw data
   * @return value_type* Pointer to data
   */
  value_type* data() noexcept { return data_; }
  /**
   * @brief Get pointer to raw data (const)
   * @return const value_type* Pointer to data
   */
  const value_type* data() const noexcept { return data_; }
  /**
   * @brief Get number of rows
   * @return size_t Number of rows
   */
  size_t rows() const noexcept { return rows_; }
  /**
   * @brief Get number of columns
   * @return size_t Number of columns
   */
  size_t cols() const noexcept { return cols_; }

  /**
   * @brief Set all elements to zero
   */
  void setZero() noexcept {
    for (size_t i = 0; i < rows_ * cols_; ++i)
      data_[i] = value_type(0);
  }

  /**
   * @brief Resize the matrix
   * @param rows New number of rows
   * @param cols New number of columns
   */
  void resize(size_t rows, size_t cols) {
    if (rows == rows_ && cols == cols_) return;
    delete[] data_;
    rows_ = rows;
    cols_ = cols;
    data_ = new value_type[rows * cols]();
  }

private:
  size_t rows_ = 0;
  size_t cols_ = 0;
  value_type* data_ = nullptr;
};

static inline CMatrix* get_matrix_ptr(C::GauXCMatrix matrix) noexcept {
  return static_cast<CMatrix*>(matrix.ptr);
}


} // namespace GauXC::detail