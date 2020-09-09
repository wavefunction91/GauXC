#include "mklsycl_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>

namespace GauXC {
namespace sycl  {
namespace blas  {

    template <typename T>
    void increment_kernel( const T* X, T* Y , sycl::nd_item<3> item_ct) {
        const auto tid = item_ct.get_group(2);
        if( tid < 1 ) (*Y) += (*X);
    }
    template <typename T> void increment(const T *X, T *Y, sycl::queue *stream) {
        stream->submit([&](sycl::handler &cgh) {
                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, 1), sycl::range<3>(1, 1, 1)),
                    [=](sycl::nd_item<3> item_ct) {
                        increment_kernel(X, Y, item_ct);
                    });
            });
    }

    template <>
    void dot(sycl::queue *handle, int N, const double *X, int INCX, const double *Y,
             int INCY, double *RES) {

        auto stat = ONEAPI::mkl::blas::dot(*handle, N, X, INCX, Y, INCY, RES);
        handle->wait();

        // if (sycl::get_pointer_type(RES, handle->get_context()) !=
        //     sycl::usm::alloc::device &&
        //     sycl::get_pointer_type(RES, handle->get_context()) !=
        //     sycl::usm::alloc::shared)
    }

    template <typename T>
    void gdot(sycl::queue *handle, int N, const T *X, int INCX, const T *Y,
              int INCY, T *SCR, T *RES) {

        dot( handle, N, X, INCX, Y, INCY, SCR );
        increment( SCR, RES, handle );
    }
    template
    void gdot( cublasHandle_t handle,
               int            N,
               const double*  X,
               int            INCX,
               const double*  Y,
               int            INCY,
               double*        SCR,
               double*        RES );


    template <typename T>
    void hadamard_product_kernel( int      M,
                                  int      N,
                                  const T* A,
                                  int      LDA,
                                  T*       B,
                                  int      LDB ,
                                  sycl::nd_item<3> item_ct) {

        const int tid_x = item_ct.get_group(2) * item_ct.get_local_range().get(2) +
            item_ct.get_local_id(2);
        const int tid_y = item_ct.get_group(1) * item_ct.get_local_range().get(1) +
            item_ct.get_local_id(1);

        if( tid_x < M and tid_y < N ) {
            B[ tid_x + tid_y*LDB ] *= A[ tid_x + tid_y*LDA ];
        }
    }

    template <typename T>
    void hadamard_product(sycl::queue *handle, int M, int N, const T *A, int LDA,
                          T *B, int LDB) {

        sycl::range<3> threads(32, 32, 1);
        sycl::range<3> blocks(util::div_ceil(M, threads[0]),
                                  util::div_ceil(N, threads[1]), 1);

        handle->submit([&](sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    sycl::nd_range<3>(
                        sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](sycl::nd_item<3> item_ct) {
                        hadamard_product_kernel(M, N, A, LDA, B, LDB, item_ct);
                    });
            });
    }
    template
    void hadamard_product( sycl::queue    *handle,
                           int            M,
                           int            N,
                           const double*  A,
                           int            LDA,
                           double*        B,
                           int            LDB );

}
}
}
