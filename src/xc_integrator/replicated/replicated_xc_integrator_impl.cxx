/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
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
                value_type* EXC, const IntegratorSettingsXC& ks_settings ) {

    eval_exc_vxc_(m,n,P,ldp,VXC,ldvxc,EXC,ks_settings);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings) {

    eval_exc_vxc_(m,n,Ps,ldps,
                      Pz,ldpz,
                      VXCs,ldvxcs,
                      VXCz,ldvxcz,EXC, ks_settings);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      const value_type* Py,
                      int64_t ldpy,
                      const value_type* Px,
                      int64_t ldpx,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* VXCy, int64_t ldvxcy,
                      value_type* VXCx, int64_t ldvxcx,
                      value_type* EXC,  const IntegratorSettingsXC& ks_settings ) {

    eval_exc_vxc_(m,n,Ps,ldps,
                      Pz,ldpz,
                      Py,ldpy,
                      Px,ldpx,
                      VXCs,ldvxcs,
                      VXCz,ldvxcz,
                      VXCy,ldvxcy,
                      VXCx,ldvxcx,EXC, ks_settings);

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_neo_exc_vxc( int64_t m1, int64_t n1, int64_t m2, int64_t n2, 
                    const value_type* P1s, int64_t ldp1s,
                    const value_type* P2s, int64_t ldp2s,
                    const value_type* P2z, int64_t ldp2z,
                    value_type* VXC1s, int64_t ldvxc1s,
                    value_type* VXC2s, int64_t ldvxc2s,
                    value_type* VXC2z, int64_t ldvxc2z,
                    value_type* EXC1,  value_type* EXC2){

    eval_neo_exc_vxc_(m1,n1,m2,n2, 
                      P1s,ldp1s,
                      P2s,ldp2s,
                      P2z,ldp2z,
                      VXC1s,ldvxc1s,
                      VXC2s,ldvxc2s,
                      VXC2z,ldvxc2z,
                      EXC1, EXC2); 

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_neo_exc_vxc( int64_t m1, int64_t n1, int64_t m2, int64_t n2, 
                    const value_type* P1s, int64_t ldp1s,
                    const value_type* P1z, int64_t ldp1z,
                    const value_type* P2s, int64_t ldp2s,
                    const value_type* P2z, int64_t ldp2z,
                    value_type* VXC1s, int64_t ldvxc1s,
                    value_type* VXC1z, int64_t ldvxc1z,
                    value_type* VXC2s, int64_t ldvxc2s,
                    value_type* VXC2z, int64_t ldvxc2z,
                    value_type* EXC1,  value_type* EXC2 ){

    eval_neo_exc_vxc_(m1,n1,m2,n2, 
                      P1s,ldp1s,
                      P1z,ldp1z,
                      P2s,ldp2s,
                      P2z,ldp2z,
                      VXC1s,ldvxc1s,
                      VXC1z,ldvxc1z,
                      VXC2s,ldvxc2s,
                      VXC2z,ldvxc2z,
                      EXC1, EXC2); 

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
