#pragma once
#include "device/type_erased_blas_handle.hpp"

namespace GauXC {

enum class DeviceBlasOp : unsigned char {
  NoTrans,
  Trans
};

enum class DeviceBlasUplo : unsigned char {
  Upper,
  Lower
};

template <typename T>
void dot( type_erased_blas_handle handle,
          int            N,
          const T*       X,
          int            INCX,
          const T*       Y,
          int            INCY,
          T*             RES );

template <typename T>
void gdot( type_erased_blas_handle handle,
          int            N,
           const T*       X,
           int            INCX,
           const T*       Y,
           int            INCY,
           T*             SCR,
           T*             RES );


template <typename T>
void hadamard_product( type_erased_blas_handle handle,
                       int            M,
                       int            N,
                       const T*       A,
                       int            LDA,
                       T*             B,
                       int            LDB );
                       

template <typename T>
void gemm( type_erased_blas_handle handle, 
           DeviceBlasOp TA, DeviceBlasOp TB,
           int M, int N, int K, T ALPHA, 
           const T* A, int LDA, const T* B, int LDB,
           T BETA, T* C, int LDC );

template <typename T>
void syr2k( type_erased_blas_handle handle, 
            DeviceBlasUplo UPLO, DeviceBlasOp Trans,
            int M, int K, T ALPHA, 
            const T* A, int LDA, const T* B, int LDB,
            T BETA, T* C, int LDC );
}

