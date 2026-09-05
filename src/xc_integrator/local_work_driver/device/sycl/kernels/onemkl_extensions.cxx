/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy).
 *
 * (c) 2024-2025, Microsoft Corporation
 *
 * All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/util/div_ceil.hpp>
#include "exceptions/onemkl_exception.hpp"

#include "device_specific/sycl_device_constants.hpp"
#include "device_specific/onemkl_util.hpp"
#include "device/common/device_blas.hpp"

// Only the BLAS domain is needed. The umbrella <oneapi/mkl.hpp> also pulls in
// the C LAPACK prototypes, whose parameter named "csl" collides with the csl
// macro from buffer_adaptor.hpp
#include <oneapi/mkl/blas.hpp>

namespace GauXC {

oneapi::mkl::transpose device_op_to_onemkl( DeviceBlasOp op ) {
  switch( op ) {
    case DeviceBlasOp::NoTrans: return oneapi::mkl::transpose::nontrans;
    case DeviceBlasOp::Trans:   return oneapi::mkl::transpose::trans;
    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported DeviceBlasOp");
      return oneapi::mkl::transpose::nontrans;
  }
}

oneapi::mkl::uplo device_uplo_to_onemkl( DeviceBlasUplo uplo ) {
  switch(uplo) {
    case DeviceBlasUplo::Upper: return oneapi::mkl::uplo::upper;
    case DeviceBlasUplo::Lower: return oneapi::mkl::uplo::lower;
    default:
      GAUXC_GENERIC_EXCEPTION("Unsupported DeviceBlasUplo");
      return oneapi::mkl::uplo::lower;
  }
}

template <typename T>
void increment( const T* X, T* Y, ::sycl::queue& stream ) {
  GAUXC_SYCL_ERROR( "Increment Launch Failed",
    stream.single_task( [=](){ (*Y) += (*X); } )
  );
}

template <typename T>
void increment( device_blas_handle generic_handle, const T* X, T* Y, int N) {
  const int threads = sycl::warp_size * sycl::max_warps_per_thread_block;
  const int blocks  = util::div_ceil( N, threads );
  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);

  GAUXC_SYCL_ERROR( "Increment Vec Launch Failed",
    stream.parallel_for( ::sycl::nd_range<1>( blocks*threads, threads ),
      [=](::sycl::nd_item<1> it) {
        const auto tid = it.get_global_id(0);
        if( tid < (size_t)N ) Y[tid] += X[tid];
      })
  );
}

template
  void increment( device_blas_handle generic_handle, const double* X, double* Y, int N );

template <>
void dot( device_blas_handle generic_handle,
          int            N,
          const double*  X,
          int            INCX,
          const double*  Y,
          int            INCY,
          double*        RES ) {

  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);

  GAUXC_ONEMKL_ERROR("ONEMKL DDOT FAILED",
    oneapi::mkl::blas::column_major::dot( stream, N, X, INCX, Y, INCY, RES ) );

}

template <typename T>
void gdot( device_blas_handle generic_handle,
           int       N,
           const T*  X,
           int       INCX,
           const T*  Y,
           int       INCY,
           T*        SCR,
           T*        RES ) {


  dot( generic_handle, N, X, INCX, Y, INCY, SCR );
  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);
  increment( SCR, RES, stream );

}

template 
void gdot( device_blas_handle generic_handle,
           int            N,
           const double*  X,
           int            INCX,
           const double*  Y,
           int            INCY,
           double*        SCR,
           double*        RES );










template <typename T>
void hadamard_product( device_blas_handle generic_handle,
                       int            M,
                       int            N,
                       const T*       A,
                       int            LDA,
                       T*             B,
                       int            LDB ) {


  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);

  // Mirrors the CUDA launch: M is carried on the fast dimension (a warp /
  // sub-group wide) and N on the slow one. SYCL orders nd_range dimensions
  // the other way round from CUDA, so the fast dimension is the last one.
  ::sycl::range<2> local( sycl::max_warps_per_thread_block, sycl::warp_size );
  ::sycl::range<2> global( util::div_ceil( N, local[0] ) * local[0],
                           util::div_ceil( M, local[1] ) * local[1] );

  GAUXC_SYCL_ERROR( "Hadamard Product Launch Failed",
    stream.parallel_for( ::sycl::nd_range<2>(global, local),
      [=](::sycl::nd_item<2> it) {
        const size_t tid_x = it.get_global_id(1);
        const size_t tid_y = it.get_global_id(0);
        if( tid_x < (size_t)M and tid_y < (size_t)N ) {
          B[ tid_x + tid_y*LDB ] *= A[ tid_x + tid_y*LDA ];
        }
      })
  );

}
 
template 
void hadamard_product( device_blas_handle generic_handle,
                       int            M,
                       int            N,
                       const double*  A,
                       int            LDA,
                       double*        B,
                       int            LDB ); 




template <>
void gemm( device_blas_handle generic_handle, 
           DeviceBlasOp TA, DeviceBlasOp TB,
           int M, int N, int K, double ALPHA, 
           const double* A, int LDA, const double* B, int LDB,
           double BETA, double* C, int LDC ) {


  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);

  GAUXC_ONEMKL_ERROR("ONEMKL DGEMM FAILED",
    oneapi::mkl::blas::column_major::gemm( stream, device_op_to_onemkl(TA),
      device_op_to_onemkl(TB), M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) );

}


template <>
void syr2k( device_blas_handle generic_handle, 
            DeviceBlasUplo UPLO, DeviceBlasOp Trans,
            int M, int K, double ALPHA, 
            const double* A, int LDA, const double* B, int LDB,
            double BETA, double* C, int LDC ) {

  auto& handle = generic_handle.blas_handle_as<util::onemkl_handle>();
  auto& stream = util::get_queue(handle);

  GAUXC_ONEMKL_ERROR("ONEMKL DSYR2K FAILED",
    oneapi::mkl::blas::column_major::syr2k( stream, device_uplo_to_onemkl(UPLO),
      device_op_to_onemkl(Trans), M, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) );

}


}
