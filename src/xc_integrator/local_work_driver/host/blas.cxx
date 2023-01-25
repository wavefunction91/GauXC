/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for retails
 */
#include "host/blas.hpp"
#include <type_traits>
#include <gauxc/exceptions.hpp>

extern "C" {

//void dlacpy_( const char* UPLO, const int* M, const int* N, const double* A, 
//              const int* LDA, double* B, const int* LDB );
//void slacpy_( const char* UPLO, const int* M, const int* N, const float* A, 
//              const int* LDA, float* B, const int* LDB );

void dgemm_( const char* TA, const char* TB, const int* M, const int* N, 
             const int* K, const double* ALPHA, const double* A, 
             const int* LDA, const double* B, const int* LDB, 
             const double* BETA, double* C, const int* LDC );
void sgemm_( const char* TA, const char* TB, const int* M, const int* N, 
             const int* K, const float* ALPHA, const float* A, 
             const int* LDA, const float* B, const int* LDB, 
             const float* BETA, float* C, const int* LDC );

void dsyr2k_( const char* UPLO, const char* TRANS, const int* N, const int* K, 
              const double* ALPHA, const double* A, const int* LDA, const double* B, 
              const int* LDB, const double* BETA, double* C, const int* LDC ); 
void ssyr2k_( const char* UPLO, const char* TRANS, const int* N, const int* K, 
              const float* ALPHA, const float* A, const int* LDA, const float* B, 
              const int* LDB, const float* BETA, float* C, const int* LDC ); 

double ddot_( const int* N, const double* X, const int* INCX, const double* Y, 
              const int* INCY );
float sdot_( const int* N, const float* X, const int* INCX, const float* Y, 
              const int* INCY );


void daxpy_( const int* N, const double* ALPHA, const double* A, const int* INCX, 
             double* Y, const int* INCY );
void saxpy_( const int* N, const float* ALPHA, const float* A, const int* INCX, 
             float* Y, const int* INCY );

void dscal_( const int* N, const double* ALPHA, const double* X, const int* INCX );
void sscal_( const int* N, const float* ALPHA, const float* X, const int* INCX ); 
}

namespace GauXC::blas {

template <typename T>
void lacpy( char UPLO, int M, int N, const T* A, int LDA, T* B,
            int LDB ) {

/*
  if constexpr ( std::is_same_v<T,float> )
    slacpy_( &UPLO, &M, &N, A, &LDA, B, &LDB );
  else if constexpr ( std::is_same_v<T,double> )
    dlacpy_( &UPLO, &M, &N, A, &LDA, B, &LDB );
  else GAUXC_GENERIC_EXCEPTION("LACPY NYI");
*/

  if( UPLO == 'L' ) {

    for( int j = 0; j < N; ++j )
    for( int i = j; i < M; ++i )
      B[i + j*LDB] = A[i + j*LDA];

  } else if( UPLO == 'U' ) {

    for( int j = 0; j <  N; ++j )
    for( int i = 0; i <= j; ++i )
      B[i + j*LDB] = A[i + j*LDA];

  } else {

    for( int j = 0; j < N; ++j )
    for( int i = 0; i < M; ++i )
      B[i + j*LDB] = A[i + j*LDA];

  }

}

template void lacpy( char UPLO, int M, int N, const float* A, int LDA, 
                     float* B, int LDB );
template void lacpy( char UPLO, int M, int N, const double* A, int LDA, 
                     double* B, int LDB );









template <typename T>
void gemm( char TA, char TB, int M, int N, int K, T ALPHA, 
           const T* A, int LDA, const T* B, int LDB, T BETA,
           T* C, int LDC ) {


  if constexpr ( std::is_same_v<T,float> )
    sgemm_( &TA, &TB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC );
  else if constexpr ( std::is_same_v<T,double> )
    dgemm_( &TA, &TB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC );
  else GAUXC_GENERIC_EXCEPTION("GEMM NYI");


}
template
void gemm( char floatA, char floatB, int M, int N, int K, float ALPHA, 
           const float* A, int LDA, const float* B, int LDB, float BETA,
           float* C, int LDC );
template
void gemm( char doubleA, char doubleB, int M, int N, int K, double ALPHA, 
           const double* A, int LDA, const double* B, int LDB, double BETA,
           double* C, int LDC );







template <typename T>
void syr2k( char UPLO, char TRANS, int N, int K, T ALPHA,
            const T* A, int LDA, const T* B, int LDB, T BETA, 
            T* C, int LDC ) {


  if constexpr ( std::is_same_v<T,float> )
    ssyr2k_( &UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC );
  else if constexpr ( std::is_same_v<T,double> )
    dsyr2k_( &UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC );
  else GAUXC_GENERIC_EXCEPTION("SYR2K NYI");


}

template
void syr2k( char UPLO, char floatRANS, int N, int K, float ALPHA,
            const float* A, int LDA, const float* B, int LDB, float BETA, 
            float* C, int LDC );
template
void syr2k( char UPLO, char doubleRANS, int N, int K, double ALPHA,
            const double* A, int LDA, const double* B, int LDB, double BETA, 
            double* C, int LDC );
            






template <typename T>
T dot( int N, const T* X, int INCX, const T* Y, int INCY ) {

  if constexpr ( std::is_same_v<T,float> )
    return sdot_(&N, X, &INCX, Y, &INCY);
  else if constexpr ( std::is_same_v<T,double> )
    return ddot_(&N, X, &INCX, Y, &INCY);
  else GAUXC_GENERIC_EXCEPTION("DOT NYI");

  return 0.;
}

template
float dot( int N, const float* X, int INCX, const float* Y, int INCY );
template
double dot( int N, const double* X, int INCX, const double* Y, int INCY );






template <typename T>
void axpy( int N, T ALPHA, const T* X, int INCX, T* Y, int INCY ) {

  if constexpr ( std::is_same_v<T,float> )
    saxpy_(&N, &ALPHA, X, &INCX, Y, &INCY );
  else if constexpr ( std::is_same_v<T,double> )
    daxpy_(&N, &ALPHA, X, &INCX, Y, &INCY );
  else GAUXC_GENERIC_EXCEPTION("AXPY NYI");

}

template
void axpy( int N, float ALPHA, const float* A, int INCX, float* Y, 
           int INCY );
template
void axpy( int N, double ALPHA, const double* A, int INCX, double* Y, 
           int INCY );
            





template <typename T>
void scal( int N, T ALPHA, T* X, int INCX ) {

  if constexpr ( std::is_same_v<T,float> )
    sscal_(&N, &ALPHA, X, &INCX );
  else if constexpr ( std::is_same_v<T,double> )
    dscal_(&N, &ALPHA, X, &INCX );
  else GAUXC_GENERIC_EXCEPTION("SCAL NYI");

}

template
void scal( int N, float ALPHA, float* X, int INCX ); 
template
void scal( int N, double ALPHA, double* X, int INCX );

}


