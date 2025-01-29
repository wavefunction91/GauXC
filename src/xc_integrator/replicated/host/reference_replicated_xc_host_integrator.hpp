/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/replicated/replicated_xc_host_integrator.hpp>
#include "xc_host_data.hpp"

namespace GauXC::detail {

template <typename ValueType>
class ReferenceReplicatedXCHostIntegrator : 
  public ReplicatedXCHostIntegrator<ValueType> {

  using base_type  = ReplicatedXCHostIntegrator<ValueType>;

public:

  static constexpr bool is_device = false;
  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;
  using task_container = std::vector<XCTask>;
  using task_iterator  = typename task_container::iterator;


protected:

  // Density Integration 
  void integrate_den_( int64_t m, int64_t n, const value_type* P, int64_t ldp, value_type* N_EL ) override;

  /// RKS EXC
  void eval_exc_( int64_t m, int64_t n, const value_type* P, int64_t ldp, 
                  value_type* EXC, const IntegratorSettingsXC& ks_settings ) override;

  /// UKS EXC
  void eval_exc_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                  const value_type* Pz, int64_t ldpz,
                  value_type* EXC, const IntegratorSettingsXC& ks_settings ) override;

  /// GKS EXC - also serves as the generic implementation
  void eval_exc_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                      const value_type* Pz, int64_t ldpz,
                      const value_type* Py, int64_t ldpy,
                      const value_type* Px, int64_t ldpx,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings ) override;

  /// RKS EXC/VXC
  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P, int64_t ldp, 
                      value_type* VXC, int64_t ldvxc, value_type* EXC, 
                      const IntegratorSettingsXC& ks_settings ) override;

  /// UKS EXC/VXC
  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                      const value_type* Pz, int64_t ldpz,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings ) override;

  /// GKS EXC/VXC - also serves as the generic implementation
  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps,
                      const value_type* Pz, int64_t ldpz,
                      const value_type* Py, int64_t ldpy,
                      const value_type* Px, int64_t ldpx,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* VXCy, int64_t ldvxcy,
                      value_type* VXCx, int64_t ldvxcx,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings ) override;


  /// RKS EXC Gradient
  void eval_exc_grad_( int64_t m, int64_t n, const value_type* P, int64_t ldp, 
                       value_type* EXC_GRAD ) override;
  /// UKS EXC Gradient
  void eval_exc_grad_( int64_t m, int64_t n, const value_type* Ps, int64_t ldps, 
                       const value_type* Pz, int64_t lpdz, value_type* EXC_GRAD ) override;

  /// sn-LinK
  void eval_exx_( int64_t m, int64_t n, const value_type* P,
                  int64_t ldp, value_type* K, int64_t ldk,
                  const IntegratorSettingsEXX& settings ) override;



  // Implementation details of integrate_den
  void integrate_den_local_work_( const value_type* P, int64_t ldp, 
                                   value_type *N_EL );

  // Implementation details of exc_vxc (for RKS/UKS/GKS deduced from input character)
  void exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            value_type* VXCs, int64_t ldvxcs,
                            value_type* VXCz, int64_t ldvxcz,
                            value_type* VXCy, int64_t ldvxcy,
                            value_type* VXCx, int64_t ldvxcx,
                            value_type* EXC, value_type *N_EL, const IntegratorSettingsXC& ks_settings,
                            task_iterator task_begin, task_iterator task_end );
                            
  // Implemetation details of exc_grad
  void exc_grad_local_work_( const value_type* Ps, int64_t ldps, const value_type* Pz, int64_t ldpz,
                             value_type* EXC_GRAD );

  // Implementation details of sn-LinK
  void exx_local_work_( const value_type* P, int64_t ldp, value_type* K, int64_t ldk,
    const IntegratorSettingsEXX& settings );

public:

  template <typename... Args>
  ReferenceReplicatedXCHostIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ReferenceReplicatedXCHostIntegrator() noexcept;


  template <typename... Args>
  void exc_vxc_local_work(Args&&... args) {
    exc_vxc_local_work_( std::forward<Args>(args)... );
  }


};

extern template class ReferenceReplicatedXCHostIntegrator<double>;

} // namespace GauXC::detail
