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
                      value_type* EXC, const IntegratorSettingsXC& ks_settings ) {

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
  neo_eval_exc_vxc( int64_t elec_m, int64_t elec_n, int64_t prot_m, int64_t prot_n, 
                    const value_type* elec_Ps, int64_t elec_ldps,
                    const value_type* prot_Ps, int64_t prot_ldps,
                    const value_type* prot_Pz, int64_t prot_ldpz,
                    value_type* elec_VXCs,     int64_t elec_ldvxcs,
                    value_type* prot_VXCs,     int64_t prot_ldvxcs,
                    value_type* prot_VXCz,     int64_t prot_ldvxcz,
                    value_type* elec_EXC,  value_type* prot_EXC,
                    const IntegratorSettingsXC& ks_settings ){

    neo_eval_exc_vxc_(elec_m,elec_n,prot_m,prot_n, 
                      elec_Ps,  elec_ldps,
                      prot_Ps,  prot_ldps,
                      prot_Pz,  prot_ldpz,
                      elec_VXCs,elec_ldvxcs,
                      prot_VXCs,prot_ldvxcs,
                      prot_VXCz,prot_ldvxcz,
                      elec_EXC, prot_EXC, ks_settings); 

}

template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  neo_eval_exc_vxc( int64_t elec_m, int64_t elec_n, int64_t prot_m, int64_t prot_n, 
                    const value_type* elec_Ps, int64_t elec_ldps,
                    const value_type* elec_Pz, int64_t elec_ldpz,
                    const value_type* prot_Ps, int64_t prot_ldps,
                    const value_type* prot_Pz, int64_t prot_ldpz,
                    value_type* elec_VXCs,     int64_t elec_ldvxcs,
                    value_type* elec_VXCz,     int64_t elec_ldvxcz,
                    value_type* prot_VXCs,     int64_t prot_ldvxcs,
                    value_type* prot_VXCz,     int64_t prot_ldvxcz,
                    value_type* elec_EXC,  value_type* prot_EXC,
                    const IntegratorSettingsXC& ks_settings ){

    neo_eval_exc_vxc_(elec_m,elec_n,prot_m,prot_n, 
                      elec_Ps,  elec_ldps,
                      elec_Pz,  elec_ldpz,
                      prot_Ps,  prot_ldps,
                      prot_Pz,  prot_ldpz,
                      elec_VXCs,elec_ldvxcs,
                      elec_VXCz,elec_ldvxcz,
                      prot_VXCs,prot_ldvxcs,
                      prot_VXCz,prot_ldvxcz,
                      elec_EXC, prot_EXC, ks_settings);  

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
