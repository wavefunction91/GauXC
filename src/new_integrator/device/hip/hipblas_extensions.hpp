#pragma once
#include <hipblas.h>

namespace GauXC {
namespace hip  {
namespace blas  {

template <typename T>
void dot( hipblasHandle_t handle,
          int            N,
          const T*       X,
          int            INCX,
          const T*       Y,
          int            INCY,
          T*             RES );

template <typename T>
void gdot( hipblasHandle_t handle,
          int            N,
           const T*       X,
           int            INCX,
           const T*       Y,
           int            INCY,
           T*             SCR,
           T*             RES );


template <typename T>
void hadamard_product( hipblasHandle_t handle,
                       int            M,
                       int            N,
                       const T*       A,
                       int            LDA,
                       T*             B,
                       int            LDB );
                       

template <typename T>
void gemm( hipblasHandle_t handle, 
           hipblasOperation_t TA, hipblasOperation_t TB,
           int M, int N, int K, T ALPHA, 
           const T* A, int LDA, const T* B, int LDB,
           T BETA, T* C, int LDC );

template <typename T>
void syr2k( hipblasHandle_t handle, 
            hipblasFillMode_t UPLO, hipblasOperation_t Trans,
            int M, int K, T ALPHA, 
            const T* A, int LDA, const T* B, int LDB,
            T BETA, T* C, int LDC );
}
}
}
