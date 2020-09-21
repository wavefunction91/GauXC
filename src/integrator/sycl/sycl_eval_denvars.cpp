#include "sycl_eval_denvars.hpp"
#include "sycl_extensions.hpp"
#include <gauxc/util/div_ceil.hpp>
#include <gauxc/exceptions/sycl_exception.hpp>

template <typename T>
using relaxed_atomic_ref = 
  cl::sycl::intel::atomic_ref< T, cl::sycl::intel::memory_order::relaxed,
                                  cl::sycl::intel::memory_scope::device,
                                  cl::sycl::access::address_space::global_space>;

namespace GauXC      {
namespace integrator {
namespace sycl       {

    template <typename T>
    void eval_uvars_lda_kernel( size_t           ntasks,
                                XCTaskDevice<T>* tasks_device ,
                                cl::sycl::nd_item<3> item_ct) {

        const size_t batch_idx = item_ct.get_group(0);
        if( batch_idx >= ntasks ) return;

        auto& task = tasks_device[ batch_idx ];

        const auto npts = task.npts;
        const auto nbf  = task.nbe;

        auto* den_eval_device   = task.den;

        const auto* basis_eval_device = task.bf;

        const auto* den_basis_prod_device = task.zmat;

        const int tid_x = item_ct.get_group(2) * item_ct.get_local_range().get(2) +
            item_ct.get_local_id(2);
        const int tid_y = item_ct.get_group(1) * item_ct.get_local_range().get(1) +
            item_ct.get_local_id(1);

        T den_reg = 0.;

        if( tid_x < nbf && tid_y < npts ) {
            const T* bf_col = basis_eval_device     + tid_y*nbf;
            const T* db_col = den_basis_prod_device + tid_y*nbf;

            den_reg = bf_col[ tid_x ]   * db_col[ tid_x ];
        }

        // Warp blocks are stored col major
        den_reg = 2 * warpReduceSum(den_reg, item_ct);

        if (item_ct.get_local_id(2) == 0 && tid_y < npts) {
            // atomic_fetch_add() to be depricated 2020 and replaced with
            // C++20 styled atomic_ref<>::fetch_add()
            // double atomicAdd(double* address, double val); // from CUDA
            // https://github.com/intel/llvm/blob/d4948b01ffbe8134ad8ebc5671e828ddbb3bd672/sycl/include/CL/sycl/atomic.hpp
            //cl::sycl::atomic_fetch_add(den_eval_device + tid_y, den_reg);
            relaxed_atomic_ref<T>( den_eval_device[tid_y] ).fetch_add( den_reg );  
        }
    }

    template <typename T>
    void eval_uvars_gga_kernel( size_t           ntasks,
                                XCTaskDevice<T>* tasks_device ,
                                cl::sycl::nd_item<3> item_ct) {

        const size_t batch_idx = item_ct.get_group(0);
        if( batch_idx >= ntasks ) return;

        auto& task = tasks_device[ batch_idx ];

        const auto npts = task.npts;
        const auto nbf  = task.nbe;

        auto* den_eval_device   = task.den;
        auto* den_x_eval_device = task.ddenx;
        auto* den_y_eval_device = task.ddeny;
        auto* den_z_eval_device = task.ddenz;

        const auto* basis_eval_device = task.bf;
        const auto* dbasis_x_eval_device = task.dbfx;
        const auto* dbasis_y_eval_device = task.dbfy;
        const auto* dbasis_z_eval_device = task.dbfz;

        const auto* den_basis_prod_device = task.zmat;

        const int tid_x = item_ct.get_group(2) * item_ct.get_local_range().get(2) +
                    item_ct.get_local_id(2);
        const int tid_y = item_ct.get_group(1) * item_ct.get_local_range().get(1) +
                    item_ct.get_local_id(1);

        T den_reg = 0.;
        T dx_reg  = 0.;
        T dy_reg  = 0.;
        T dz_reg  = 0.;

        if( tid_x < nbf && tid_y < npts ) {

            const T* bf_col   = basis_eval_device     + tid_y*nbf;
            const T* bf_x_col = dbasis_x_eval_device  + tid_y*nbf;
            const T* bf_y_col = dbasis_y_eval_device  + tid_y*nbf;
            const T* bf_z_col = dbasis_z_eval_device  + tid_y*nbf;
            const T* db_col   = den_basis_prod_device + tid_y*nbf;

            den_reg = bf_col[ tid_x ]   * db_col[ tid_x ];
            dx_reg  = bf_x_col[ tid_x ] * db_col[ tid_x ];
            dy_reg  = bf_y_col[ tid_x ] * db_col[ tid_x ];
            dz_reg  = bf_z_col[ tid_x ] * db_col[ tid_x ];

        }

        // Warp blocks are stored col major
        den_reg = 2 * warpReduceSum(den_reg, item_ct);
        dx_reg = 4 * warpReduceSum(dx_reg, item_ct);
        dy_reg = 4 * warpReduceSum(dy_reg, item_ct);
        dz_reg = 4 * warpReduceSum(dz_reg, item_ct);

        if (item_ct.get_local_id(2) == 0 && tid_y < npts) {
            //cl::sycl::atomic_fetch_add(den_eval_device + tid_y, den_reg);
            //cl::sycl::atomic_fetch_add(den_x_eval_device + tid_y, dx_reg);
            //cl::sycl::atomic_fetch_add(den_y_eval_device + tid_y, dy_reg);
            //cl::sycl::atomic_fetch_add(den_z_eval_device + tid_y, dz_reg);
            relaxed_atomic_ref<T>( den_eval_device[tid_y] ).  fetch_add( den_reg );  
            relaxed_atomic_ref<T>( den_x_eval_device[tid_y] ).fetch_add( dx_reg );  
            relaxed_atomic_ref<T>( den_y_eval_device[tid_y] ).fetch_add( dy_reg );  
            relaxed_atomic_ref<T>( den_z_eval_device[tid_y] ).fetch_add( dz_reg );  
        }
    }


    template <typename T>
    void eval_vvars_gga_kernel(size_t   npts,
                               const T* den_x_eval_device,
                               const T* den_y_eval_device,
                               const T* den_z_eval_device,
                               T* gamma_eval_device,
                               cl::sycl::nd_item<3> item_ct) {

        const size_t tid = item_ct.get_local_id(2) + item_ct.get_group(2) * item_ct.get_local_range().get(2);

        if( tid < npts ) {
            const T dx = den_x_eval_device[ tid ];
            const T dy = den_y_eval_device[ tid ];
            const T dz = den_z_eval_device[ tid ];

            gamma_eval_device[tid] = dx*dx + dy*dy + dz*dz;
        }
    }

    template <typename T>
    void eval_uvars_lda_device(size_t ntasks, size_t max_nbf, size_t max_npts,
                               XCTaskDevice<T> *tasks_device, cl::sycl::queue *queue) {

        cl::sycl::range<3> threads(32, 32, 1);
        cl::sycl::range<3> blocks(util::div_ceil(max_nbf, 32),
                                  util::div_ceil(max_npts, 32),
                                  ntasks);

        GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(cl::sycl::nd_range<3>(cl::sycl::range<3>(global_range.get(2),
                                                                          global_range.get(1),
                                                                          global_range.get(0)),
                                                       cl::sycl::range<3>(threads.get(2),
                                                                          threads.get(1),
                                                                          threads.get(0))),

                    [=](cl::sycl::nd_item<3> item_ct) {
                        eval_uvars_lda_kernel(ntasks, tasks_device, item_ct);
                    });
                }) );
    }

    template <typename T>
    void eval_uvars_gga_device(size_t ntasks, size_t max_nbf, size_t max_npts,
                               XCTaskDevice<T> *tasks_device, cl::sycl::queue *queue) {

        cl::sycl::range<3> threads(32, 32, 1);
        cl::sycl::range<3> blocks(util::div_ceil(max_nbf, 32),
                                  util::div_ceil(max_npts, 32),
                                  ntasks);

        GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(cl::sycl::nd_range<3>(cl::sycl::range<3>(global_range.get(2),
                                                                          global_range.get(1),
                                                                          global_range.get(0)),
                                                       cl::sycl::range<3>(threads.get(2),
                                                                          threads.get(1),
                                                                          threads.get(0))),

                    [=](cl::sycl::nd_item<3> item_ct) {
                        eval_uvars_gga_kernel(ntasks, tasks_device, item_ct);
                    });
                }) );
    }

    template <typename T>
    void eval_vvars_gga_device(size_t npts,
                               const T *den_x_device,
                               const T *den_y_device,
                               const T *den_z_device,
                               T *gamma_device,
                               cl::sycl::queue *queue) {

        cl::sycl::range<3> threads(1024, 1, 1);
        cl::sycl::range<3> blocks(util::div_ceil(npts, 1024), 1, 1);

        GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
                auto global_range = blocks * threads;

                cgh.parallel_for(cl::sycl::nd_range<3>(cl::sycl::range<3>(global_range.get(2),
                                                                          global_range.get(1),
                                                                          global_range.get(0)),
                                                       cl::sycl::range<3>(threads.get(2),
                                                                          threads.get(1),
                                                                          threads.get(0))),

                    [=](cl::sycl::nd_item<3> item_ct) {
                        eval_vvars_gga_kernel(npts, den_x_device, den_y_device, den_z_device,
                                              gamma_device, item_ct);
                    });
                }) );
    }





template
void eval_uvars_lda_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            cl::sycl::queue*          stream );

template
void eval_uvars_gga_device( size_t                ntasks,
                            size_t                max_nbf,
                            size_t                max_npts,
                            XCTaskDevice<double>* tasks_device,
                            cl::sycl::queue*          stream );

template
void eval_vvars_gga_device( size_t            npts,
                            const double*     den_x_device,
                            const double*     den_y_device,
                            const double*     den_z_device,
                                  double*     gamma_device,
                            cl::sycl::queue*      stream );
}
}
}
