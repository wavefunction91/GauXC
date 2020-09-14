#include "mklsycl_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/exceptions/sycl_exception.hpp>

namespace GauXC {
namespace sycl  {
namespace blas  {

    template <typename T>
    void increment_kernel( const T* X, T* Y , cl::sycl::nd_item<3> item_ct) {
        const auto tid = item_ct.get_group(2);
        if( tid < 1 ) (*Y) += (*X);
    }
    template <typename T> void increment(const T *X, T *Y, cl::sycl::queue *queue) {
        GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
                cgh.parallel_for(
                    cl::sycl::nd_range<3>(cl::sycl::range<3>(1, 1, 1), cl::sycl::range<3>(1, 1, 1)),
                    [=](cl::sycl::nd_item<3> item_ct) {
                        increment_kernel(X, Y, item_ct);
                    });
                }) );
    }

    template <>
    void dot(cl::sycl::queue *syclQue, int N, const double *X, int INCX, const double *Y,
             int INCY, double *RES) {

        // NOTE abb: might need to change `oneapi` to `ONEAPI` in future
        GAUXC_SYCL_ERROR( auto stat = oneapi::mkl::blas::dot(*syclQue, N, X, INCX, Y, INCY, RES) );
        syclQue->wait();  // abb: get rid of wait() after build, test

        // if (cl::sycl::get_pointer_type(RES, syclQue->get_context()) !=
        //     cl::sycl::usm::alloc::device &&
        //     cl::sycl::get_pointer_type(RES, syclQue->get_context()) !=
        //     cl::sycl::usm::alloc::shared)
    }

    template <typename T>
    void gdot(cl::sycl::queue *syclQue, int N, const T *X, int INCX, const T *Y,
              int INCY, T *SCR, T *RES) {

        dot( syclQue, N, X, INCX, Y, INCY, SCR );
        increment( SCR, RES, syclQue );
    }
    template
    void gdot( cublasSyclQue_t syclQue,
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
                                  cl::sycl::nd_item<3> item_ct) {

        const int tid_x = item_ct.get_group(2) * item_ct.get_local_range().get(2) +
            item_ct.get_local_id(2);
        const int tid_y = item_ct.get_group(1) * item_ct.get_local_range().get(1) +
            item_ct.get_local_id(1);

        if( tid_x < M and tid_y < N ) {
            B[ tid_x + tid_y*LDB ] *= A[ tid_x + tid_y*LDA ];
        }
    }

    template <typename T>
    void hadamard_product(cl::sycl::queue *syclQue, int M, int N, const T *A, int LDA,
                          T *B, int LDB) {

        cl::sycl::range<3> threads(32, 32, 1);
        cl::sycl::range<3> blocks(util::div_ceil(M, threads[0]),
                                  util::div_ceil(N, threads[1]), 1);

        GAUXC_SYCL_ERROR( syclQue->submit([&](cl::sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(
                    cl::sycl::nd_range<3>(
                        cl::sycl::range<3>(global_range.get(2), global_range.get(1),
                                           global_range.get(0)),
                        cl::sycl::range<3>(threads.get(2), threads.get(1), threads.get(0))),
                    [=](cl::sycl::nd_item<3> item_ct) {
                        hadamard_product_kernel(M, N, A, LDA, B, LDB, item_ct);
                    });
                }) );
    }
    template
    void hadamard_product( cl::sycl::queue    *syclQue,
                           int            M,
                           int            N,
                           const double*  A,
                           int            LDA,
                           double*        B,
                           int            LDB );

}
}
}
