/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
PGASDistributedXCIntegratorImpl<ValueType>::
  PGASDistributedXCIntegratorImpl( std::shared_ptr< functional_type >   func,
                              std::shared_ptr< LoadBalancer >      lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd) :
    func_(func), load_balancer_(lb), local_work_driver_(std::move(lwd)) { }

template <typename ValueType>
PGASDistributedXCIntegratorImpl<ValueType>::
  ~PGASDistributedXCIntegratorImpl() noexcept = default;

template <typename ValueType>
void PGASDistributedXCIntegratorImpl<ValueType>::
  integrate_den( const matrix_type& P, value_type* N_EL ) {

    integrate_den_(P,N_EL);

}

template <typename ValueType>
void PGASDistributedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( const matrix_type& P, matrix_type& VXC, value_type* EXC ) {

    eval_exc_vxc_(P,VXC,EXC);

}

template <typename ValueType>
void PGASDistributedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( const matrix_type& Ps, const matrix_type& Pz,
                      matrix_type& VXCs, matrix_type& VXCz,
                      value_type* EXC ) {

    eval_exc_vxc_(Ps,Pz,VXCs,VXCz,EXC);

}

template <typename ValueType>
void PGASDistributedXCIntegratorImpl<ValueType>::
  eval_exc_grad( const matrix_type& P, value_type* EXC_GRAD ) {

    eval_exc_grad_(P,EXC_GRAD);

}

template <typename ValueType>
void PGASDistributedXCIntegratorImpl<ValueType>::
  eval_exx( const matrix_type& P, matrix_type& K, 
            const IntegratorSettingsEXX& settings ) {

    eval_exx_(P,K,settings);

}

template class PGASDistributedXCIntegratorImpl<double>;

}
}
