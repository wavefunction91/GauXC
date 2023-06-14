/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/xc_integrator/replicated/replicated_xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type >   func,
                              std::shared_ptr< LoadBalancer >      lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd,
                              std::shared_ptr< ReductionDriver >   rd) :
    func_(func), load_balancer_(lb), local_work_driver_(std::move(lwd)),
    reduction_driver_(rd){ }

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ~ReplicatedXCIntegratorImpl() noexcept = default;

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  integrate_den( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* N_EL ) {

    integrate_den_(m,n,P,ldp,N_EL);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( int64_t m, int64_t n, const value_type* P,
                int64_t ldp, value_type* VXC, int64_t ldvxc,
                value_type* EXC ) {

    eval_exc_vxc_(m,n,P,ldp,VXC,ldvxc,EXC);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exc_grad( int64_t m, int64_t n, const value_type* P,
                int64_t ldp, value_type* EXC_GRAD ) {

    eval_exc_grad_(m,n,P,ldp,EXC_GRAD);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exx( int64_t m, int64_t n, const value_type* P,
            int64_t ldp, value_type* K, int64_t ldk,
            const IntegratorSettingsEXX& settings ) {

    eval_exx_(m,n,P,ldp,K,ldk,settings);

}

template class ReplicatedXCIntegratorImpl<double>;

}
}
