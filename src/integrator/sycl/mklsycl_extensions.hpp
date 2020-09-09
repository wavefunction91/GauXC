#pragma once

#include <SYCL/sycl.hpp>
#include "mkl_blas_sycl.hpp"
#include "mkl_lapack_sycl.hpp"

namespace GauXC {
namespace sycl  {
namespace blas  {

template <typename T>
void dot(sycl::queue *handle, int N, const T *X, int INCX, const T *Y, int INCY,
         T *RES);

template <typename T>
void gdot(sycl::queue *handle, int N, const T *X, int INCX, const T *Y,
          int INCY, T *SCR, T *RES);

template <typename T>
void hadamard_product(sycl::queue *handle, int M, int N, const T *A, int LDA,
                      T *B, int LDB);

}
}
}
