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
                       

template <typename T>
void gemm( cublasHandle_t handle, 
           cublasOperation_t TA, cublasOperation_t TB,
           int M, int N, int K, T ALPHA, 
           const T* A, int LDA, const T* B, int LDB,
           T BETA, T* C, int LDC );

template <typename T>
void syr2k( cublasHandle_t handle, 
            cublasFillMode_t UPLO, cublasOperation_t Trans,
            int M, int K, T ALPHA, 
            const T* A, int LDA, const T* B, int LDB,
            T BETA, T* C, int LDC );
}
}
}
