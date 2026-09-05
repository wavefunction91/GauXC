#pragma once
#include <cstdint>

namespace GauXC::blas {

template <typename T>
void lacpy( char UPLO, int M, int N, const T* A, int LDA, T* B,
            int LDB );

template <typename T>
void gemm( char TA, char TB, int M, int N, int K, T ALPHA, 
           const T* A, int LDA, const T* B, int LDB, T BETA,
           T* C, int LDC );

template <typename T>
void syr2k( char UPLO, char TRANS, int N, int K, T ALPHA,
            const T* A, int LDA, const T* B, int LDB, T BETA, 
            T* C, int LDC ); 
            

template <typename T>
T dot( int N, const T* X, int INCX, const T* Y, int INCY );

template <typename T>
void axpy( int N, T ALPHA, const T* X, int INCX, T* Y, int INCY );
            
template <typename T>
void scal( int N, T ALPHA,  T* X, int INCX );

}
