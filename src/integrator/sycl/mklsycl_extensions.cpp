#include "mklsycl_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/exceptions/sycl_exception.hpp>

namespace GauXC {
namespace sycl  {
namespace blas  {

    template <typename T> void increment(const T *X, T *Y, cl::sycl::queue *queue) {
        GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
           cgh.single_task([=] () {
               (*Y) += (*X);
           });
        }) );
    }

    template <>
    void dot(cl::sycl::queue *syclQue, int N, const double *X, int INCX, const double *Y,
             int INCY, double *RES) {

        // abb: might need to change `oneapi` to `ONEAPI` in future
        GAUXC_SYCL_ERROR( auto stat = oneapi::mkl::blas::dot(*syclQue, N, X, INCX, Y, INCY, RES) );
    }

    template <typename T>
    void gdot(cl::sycl::queue *syclQue, int N, const T *X, int INCX, const T *Y,
              int INCY, T *SCR, T *RES) {

        dot( syclQue, N, X, INCX, Y, INCY, SCR );
        increment( SCR, RES, syclQue );
    }
    template
    void gdot( cl::sycl::queue *syclQue,
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
                                  cl::sycl::item<2>& item_ct) {

        const int tid_x = item_ct.get_id(1);
        const int tid_y = item_ct.get_id(0);

        if( tid_x < M and tid_y < N ) {
            B[ tid_x + tid_y*LDB ] *= A[ tid_x + tid_y*LDA ];
        }
    }

    template <typename T>
    void hadamard_product(cl::sycl::queue *syclQue, int M, int N, const T *A, int LDA,
                          T *B, int LDB) {

        GAUXC_SYCL_ERROR( syclQue->submit([&](cl::sycl::handler &cgh) {
                    cgh.parallel_for( cl::sycl::range<2>{N, M}, [=](cl::sycl::item<2> item_ct) {
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
