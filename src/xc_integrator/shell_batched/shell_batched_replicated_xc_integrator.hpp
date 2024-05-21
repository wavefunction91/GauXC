/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/gauxc_config.hpp>
#include "shell_batched_xc_integrator.hpp"
#ifdef GAUXC_HAS_DEVICE
#include "device/xc_device_data.hpp"
#endif

namespace GauXC {
namespace detail {

template <typename BaseIntegratorType, typename IncoreIntegratorType>
class ShellBatchedReplicatedXCIntegrator : 
  public BaseIntegratorType,
  public ShellBatchedXCIntegratorBase {

  using base_type  = BaseIntegratorType;

public:

  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

  using host_task_container = std::vector<XCTask>;
  using host_task_iterator  = typename host_task_container::iterator;

protected:

#ifdef GAUXC_HAS_DEVICE
  std::unique_ptr<XCDeviceData> device_data_ptr_;
#endif

  using incore_integrator_type = IncoreIntegratorType;
  using incore_task_data = ShellBatchedXCIntegratorBase::incore_task_data;

  // Density Integration 
  void integrate_den_( int64_t m, int64_t n, const value_type* P, int64_t ldp, value_type* N_EL ) override;

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
  void eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                       int64_t ldp, value_type* EXC_GRAD ) override;

  /// sn-LinK
  void eval_exx_( int64_t m, int64_t n, const value_type* P,
                  int64_t ldp, value_type* K, int64_t ldk,
                  const IntegratorSettingsEXX& settings ) override;




  // Implementation details of exc_vxc (for RKS/UKS/GKS deduced from input character)
  void exc_vxc_local_work_( const basis_type& basis, const value_type* Ps, int64_t ldps,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Py, int64_t ldpy,
                            const value_type* Px, int64_t ldpx,
                            value_type* VXCs, int64_t ldvxcs,
                            value_type* VXCz, int64_t ldvxcz,
                            value_type* VXCy, int64_t ldvxcy,
                            value_type* VXCx, int64_t ldvxcx,
                            value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end, incore_integrator_type& incore_integrator
                             );


  void execute_task_batch( incore_task_data& task, const basis_type& basis, const Molecule& mol, 
                           const value_type* Ps, int64_t ldps,
                           const value_type* Pz, int64_t ldpz,
                           const value_type* Py, int64_t ldpy,
                           const value_type* Px, int64_t ldpx,
                           value_type* VXCs, int64_t ldvxcs,
                           value_type* VXCz, int64_t ldvxcz,
                           value_type* VXCy, int64_t ldvxcy,
                           value_type* VXCx, int64_t ldvxcx,
                           value_type* EXC, value_type* N_EL, incore_integrator_type& incore_integrator);
public:

  template <typename... Args>
  ShellBatchedReplicatedXCIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ShellBatchedReplicatedXCIntegrator() noexcept = default;

};

}
}
