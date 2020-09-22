#include "sycl_inc_potential.hpp"
#include <gauxc/exceptions/sycl_exception.hpp>

namespace GauXC      {
namespace integrator {
namespace sycl       {

template <typename T>
void inc_by_submat_combined_kernel( size_t           ntasks,
                                    XCTaskDevice<T>* device_tasks,
                                    T*               A,
                                    size_t           LDA ,
                                    cl::sycl::nd_item<3> item_ct) {

    const size_t batch_id = item_ct.get_group(0);

    if( batch_id < ntasks ) {

        auto& task = device_tasks[ batch_id ];

        const auto  ncut              = task.ncut;
        const auto* submat_cut_device = task.submat_cut;
        const auto  LDAS              = task.nbe;
        auto* ASmall_device           = task.nbe_scr;

        //if( LDAS == LDAB ) return;

        const int tid_x =
            item_ct.get_local_range().get(2) * item_ct.get_group(2) +
            item_ct.get_local_id(2);
        const int tid_y =
            item_ct.get_local_range().get(1) * item_ct.get_group(1) +
            item_ct.get_local_id(1);

        int64_t i(0);
        for( size_t i_cut = 0; i_cut < ncut; ++i_cut ) {
            const int64_t i_cut_first  = submat_cut_device[ 2*i_cut ];
            const int64_t i_cut_second = submat_cut_device[ 2*i_cut + 1 ];
            const int64_t delta_i      = i_cut_second - i_cut_first;

            int64_t j(0);
            for( size_t j_cut = 0; j_cut < ncut; ++j_cut ) {
                const int64_t j_cut_first  = submat_cut_device[ 2*j_cut ];
                const int64_t j_cut_second = submat_cut_device[ 2*j_cut + 1 ];
                const int64_t delta_j      = j_cut_second - j_cut_first;

                auto* ASmall_begin = ASmall_device + i           + j          *LDAS;
                auto* ABig_begin   = A             + i_cut_first + j_cut_first*LDA ;

                for( int64_t J = tid_y; J < delta_j; J += item_ct.get_local_range().get(1) )
                    for( int64_t I = tid_x; I < delta_i; I += item_ct.get_local_range().get(2) ) {
                        //cl::sycl::atomic_fetch_add( ABig_begin + I + J*LDA, ASmall_begin[I+J*LDAS] );
                        auto atm = cl::sycl::intel::atomic_ref<T, cl::sycl::intel::memory_order::relaxed,
                                                               cl::sycl::intel::memory_scope::device, cl::sycl::access::address_space::global_space>( *(ABig_begin + I + J*LDA) );
                        atm.fetch_add( ASmall_begin[I+J*LDAS] );
                    }
                j += delta_j;
            }
            i += delta_i;
        }

    } // batch_id check
}

template <typename T>
void task_inc_potential(size_t ntasks, XCTaskDevice<T> *device_tasks,
                        T *V_device, size_t LDV, cl::sycl::queue *queue) {

    cl::sycl::range<3> threads(16, 16, 1), blocks(1, 1, ntasks);
    GAUXC_SYCL_ERROR( queue->submit([&](cl::sycl::handler &cgh) {
            auto global_range = blocks * threads;

            cgh.parallel_for(cl::sycl::nd_range<3>(cl::sycl::range<3>(global_range.get(2),
                                                                      global_range.get(1),
                                                                      global_range.get(0)),
                                                   cl::sycl::range<3>(threads.get(2),
                                                                      threads.get(1),
                                                                      threads.get(0))),

                [=](cl::sycl::nd_item<3> item_ct) {
                    inc_by_submat_combined_kernel(ntasks, device_tasks, V_device, LDV,
                                                  item_ct);
                });
            }) );
}
template
void task_inc_potential( size_t                ntasks,
                         XCTaskDevice<double>* device_tasks,
                         double*               V_device,
                         size_t                LDV,
                         cl::sycl::queue*      queue );

}
}
}
