#pragma once
#include <cublas_v2.h>

namespace GauXC {
namespace cuda  {
namespace blas  {

template <typename T>
void dot( cublasHandle_t handle,
          int            N,
          const T*       X,
          int            INCX,
          const T*       Y,
          int            INCY,
          T*             RES );

template <typename T>
void gdot( cublasHandle_t handle,
          int            N,
           const T*       X,
           int            INCX,
           const T*       Y,
           int            INCY,
           T*             SCR,
           T*             RES );


template <typename T>
void hadamard_product( cublasHandle_t handle,
                       int            M,
                       int            N,
                       const T*       A,
                       int            LDA,
                       T*             B,
                       int            LDB );
                       
}
}
}
