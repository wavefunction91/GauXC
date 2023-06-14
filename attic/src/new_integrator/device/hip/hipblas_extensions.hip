#include "hip/hip_runtime.h"
#include "device/hip/hipblas_extensions.hpp"
#include <gauxc/util/hipblas_util.hpp>
#include <gauxc/util/div_ceil.hpp>
#include "exceptions/hipblas_exception.hpp"

#include "device/hip/hip_device_properties.hpp"

namespace GauXC {
namespace hip  {
namespace blas  {

using namespace GauXC::hip;

template <typename T>
__global__ void increment_kernel( const T* X, T* Y ) {
  const auto tid = blockIdx.x;
  if( tid < 1 ) (*Y) += (*X);
}

template <typename T>
void increment( const T* X, T* Y, hipStream_t stream ) {
  hipLaunchKernelGGL(increment_kernel, dim3(1), dim3(1), 0, stream, X,Y);
}

template <>
void dot( hipblasHandle_t handle,
          int            N,
          const double*  X,
          int            INCX,
          const double*  Y,
          int            INCY,
          double*        RES ) {

  auto stat = hipblasDdot( handle, N, X, INCX, Y, INCY, RES );
  GAUXC_HIPBLAS_ERROR("HIPBLAS DDOT FAILED", stat );

}

template <typename T>
void gdot( hipblasHandle_t handle,
           int       N,
           const T*  X,
           int       INCX,
           const T*  Y,
           int       INCY,
           T*        SCR,
           T*        RES ) {

  dot( handle, N, X, INCX, Y, INCY, SCR );
  auto stream = util::get_stream(handle);
  increment( SCR, RES, stream );

}

template 
void gdot( hipblasHandle_t handle,
           int            N,
           const double*  X,
           int            INCX,
           const double*  Y,
           int            INCY,
           double*        SCR,
           double*        RES );










template <typename T>
void __global__ hadamard_product_kernel( int      M,
                                         int      N,
                                         const T* A,
                                         int      LDA,
                                         T*       B,
                                         int      LDB ) {

  const int tid_x = blockIdx.x * blockDim.x + threadIdx.x;
  const int tid_y = blockIdx.y * blockDim.y + threadIdx.y;

  if( tid_x < M and tid_y < N ) {
    B[ tid_x + tid_y*LDB ] *= A[ tid_x + tid_y*LDA ];
  }

}



template <typename T>
void hadamard_product( hipblasHandle_t handle,
                       int            M,
                       int            N,
                       const T*       A,
                       int            LDA,
                       T*             B,
                       int            LDB ) {

  auto stream = util::get_stream(handle);
  dim3 threads(warp_size, max_warps_per_thread_block);
  dim3 blocks( util::div_ceil( M, threads.x ),
               util::div_ceil( N, threads.y ) );

  hipLaunchKernelGGL(hadamard_product_kernel, dim3(blocks), dim3(threads), 0, stream ,  M, N, A, LDA, B, LDB );

}
 
template 
void hadamard_product( hipblasHandle_t handle,
                       int            M,
                       int            N,
                       const double*  A,
                       int            LDA,
                       double*        B,
                       int            LDB ); 




template <>
void gemm( hipblasHandle_t handle, 
           hipblasOperation_t TA, hipblasOperation_t TB,
           int M, int N, int K, double ALPHA, 
           const double* A, int LDA, const double* B, int LDB,
           double BETA, double* C, int LDC ) {

  auto stat = hipblasDgemm( handle, TA, TB, M, N, K, &ALPHA, A, LDA,
                           B, LDB, &BETA, C, LDC );
  GAUXC_HIPBLAS_ERROR("HIPBLAS DGEMM FAILED", stat);

}


template <>
void syr2k( hipblasHandle_t handle, 
            hipblasFillMode_t UPLO, hipblasOperation_t Trans,
            int M, int K, double ALPHA, 
            const double* A, int LDA, const double* B, int LDB,
            double BETA, double* C, int LDC ) {

  auto stat = hipblasDsyr2k( handle, UPLO, Trans, M, K, &ALPHA, A, LDA, B, LDB,
                           &BETA, C, LDC );
  GAUXC_HIPBLAS_ERROR("HIPBLAS DSYR2K FAILED", stat);

}

}
}
}

