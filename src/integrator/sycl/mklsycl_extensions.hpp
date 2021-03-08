#pragma once

#include <CL/sycl.hpp>
#include <oneapi/mkl.hpp>

namespace GauXC {
namespace sycl  {
namespace blas  {

template <typename T>
void dot(cl::sycl::queue *syclQue, int N, const T *X, int INCX, const T *Y, int INCY,
         T *RES);

template <typename T>
void gdot(cl::sycl::queue *syclQue, int N, const T *X, int INCX, const T *Y,
          int INCY, T *SCR, T *RES);

template <typename T>
void hadamard_product(cl::sycl::queue *syclQue, int M, int N, const T *A, int LDA,
                      T *B, int LDB);

}
}
}
