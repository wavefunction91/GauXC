#include "blas.hpp"
#include <gauxc/gauxc_config.hpp>
#include <type_traits>
#include <stdexcept>

#ifdef GAUXC_BLAS_ILP64
  using blas_int = int64_t;
#else
  using blas_int = int32_t;
#endif


extern "C" {


void dgemm_( const char* TA, const char* TB, const blas_int* M, const blas_int* N, 
             const blas_int* K, const double* ALPHA, const double* A, 
             const blas_int* LDA, const double* B, const blas_int* LDB, 
             const double* BETA, double* C, const blas_int* LDC );
void sgemm_( const char* TA, const char* TB, const blas_int* M, const blas_int* N, 
             const blas_int* K, const float* ALPHA, const float* A, 
             const blas_int* LDA, const float* B, const blas_int* LDB, 
             const float* BETA, float* C, const blas_int* LDC );

void dsyr2k_( const char* UPLO, const char* TRANS, const blas_int* N, const blas_int* K, 
              const double* ALPHA, const double* A, const blas_int* LDA, const double* B, 
              const blas_int* LDB, const double* BETA, double* C, const blas_int* LDC ); 
void ssyr2k_( const char* UPLO, const char* TRANS, const blas_int* N, const blas_int* K, 
              const float* ALPHA, const float* A, const blas_int* LDA, const float* B, 
              const blas_int* LDB, const float* BETA, float* C, const blas_int* LDC ); 

double ddot_( const blas_int* N, const double* X, const blas_int* INCX, const double* Y, 
              const blas_int* INCY );
float sdot_( const blas_int* N, const float* X, const blas_int* INCX, const float* Y, 
              const blas_int* INCY );


void daxpy_( const blas_int* N, const double* ALPHA, const double* A, const blas_int* INCX, 
             double* Y, const blas_int* INCY );
void saxpy_( const blas_int* N, const float* ALPHA, const float* A, const blas_int* INCX, 
             float* Y, const blas_int* INCY );

void dscal_( const blas_int* N, const double* ALPHA, const double* X, const blas_int* INCX );
void sscal_( const blas_int* N, const float* ALPHA, const float* X, const blas_int* INCX ); 
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
  else throw std::runtime_error("LACPY NYI");
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

  blas_int _M = M;
  blas_int _N = N;
  blas_int _K = K;
  blas_int _LDA = LDA;
  blas_int _LDB = LDB;
  blas_int _LDC = LDC;

  if constexpr ( std::is_same_v<T,float> )
    sgemm_( &TA, &TB, &_M, &_N, &_K, &ALPHA, A, &_LDA, B, &_LDB, &BETA, C, &_LDC );
  else if constexpr ( std::is_same_v<T,double> )
    dgemm_( &TA, &TB, &_M, &_N, &_K, &ALPHA, A, &_LDA, B, &_LDB, &BETA, C, &_LDC );
  else throw std::runtime_error("GEMM NYI");


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

  blas_int _N = N;
  blas_int _K = K;
  blas_int _LDA = LDA;
  blas_int _LDB = LDB;
  blas_int _LDC = LDC;

  if constexpr ( std::is_same_v<T,float> )
    ssyr2k_( &UPLO, &TRANS, &_N, &_K, &ALPHA, A, &_LDA, B, &_LDB, &BETA, C, &_LDC );
  else if constexpr ( std::is_same_v<T,double> )             
    dsyr2k_( &UPLO, &TRANS, &_N, &_K, &ALPHA, A, &_LDA, B, &_LDB, &BETA, C, &_LDC );
  else throw std::runtime_error("SYR2K NYI");


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

  blas_int _N = N;
  blas_int _INCX = INCX;
  blas_int _INCY = INCY;
 
  if constexpr ( std::is_same_v<T,float> )
    return sdot_(&_N, X, &_INCX, Y, &_INCY);
  else if constexpr ( std::is_same_v<T,double> )
    return ddot_(&_N, X, &_INCX, Y, &_INCY);
  else throw std::runtime_error("DOT NYI");

  return 0.;
}

template
float dot( int N, const float* X, int INCX, const float* Y, int INCY );
template
double dot( int N, const double* X, int INCX, const double* Y, int INCY );






template <typename T>
void axpy( int N, T ALPHA, const T* X, int INCX, T* Y, int INCY ) {

  blas_int _N = N;
  blas_int _INCX = INCX;
  blas_int _INCY = INCY;

  if constexpr ( std::is_same_v<T,float> )
    saxpy_(&_N, &ALPHA, X, &_INCX, Y, &_INCY );
  else if constexpr ( std::is_same_v<T,double> )
    daxpy_(&_N, &ALPHA, X, &_INCX, Y, &_INCY );
  else throw std::runtime_error("AXPY NYI");

}

template
void axpy( int N, float ALPHA, const float* A, int INCX, float* Y, 
           int INCY );
template
void axpy( int N, double ALPHA, const double* A, int INCX, double* Y, 
           int INCY );
            





template <typename T>
void scal( int N, T ALPHA, T* X, int INCX ) {

  blas_int _N = N;
  blas_int _INCX = INCX;

  if constexpr ( std::is_same_v<T,float> )
    sscal_(&_N, &ALPHA, X, &_INCX );
  else if constexpr ( std::is_same_v<T,double> )
    dscal_(&_N, &ALPHA, X, &_INCX );
  else throw std::runtime_error("SCAL NYI");

}

template
void scal( int N, float ALPHA, float* X, int INCX ); 
template
void scal( int N, double ALPHA, double* X, int INCX );

}


