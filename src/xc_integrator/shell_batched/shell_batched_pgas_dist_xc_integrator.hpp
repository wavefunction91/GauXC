/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/gauxc_config.hpp>
#include "shell_batched_xc_integrator.hpp"
#ifdef GAUXC_ENABLE_DEVICE
#include "device/xc_device_data.hpp"
#endif

namespace GauXC {
namespace detail {

template <typename BaseIntegratorType, typename IncoreIntegratorType>
class ShellBatchedPGASDistributedXCIntegrator : 
  public BaseIntegratorType,
  public ShellBatchedXCIntegratorBase {

  using base_type  = BaseIntegratorType;

public:

  using value_type = typename base_type::value_type;
  using matrix_type = typename base_type::matrix_type;
  using basis_type = typename base_type::basis_type;

  using host_task_container = std::vector<XCTask>;
  using host_task_iterator  = typename host_task_container::iterator;

protected:

#ifdef GAUXC_ENABLE_DEVICE
  std::unique_ptr<XCDeviceData> device_data_ptr_;
#endif

  using incore_integrator_type = IncoreIntegratorType;
  using incore_task_data = ShellBatchedXCIntegratorBase::incore_task_data;

  void integrate_den_( const matrix_type& P, value_type* N_EL ) override;
  void eval_exc_vxc_( const matrix_type& P, matrix_type& VXC, value_type* EXC ) override;
  void eval_exc_vxc_( const matrix_type& Ps, const matrix_type& Pz, 
                      matrix_type& VXCs, matrix_type& VXCz, value_type* EXC ) override;
  void eval_exc_grad_( const matrix_type& P, value_type* EXC_GRAD ) override;
  void eval_exx_( const matrix_type& P, matrix_type& K,
                  const IntegratorSettingsEXX& settings ) override;

  void exc_vxc_local_work_( const basis_type& basis, const matrix_type& P, 
                            matrix_type& VXC, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            incore_integrator_type& incore_integrator );

  void execute_task_batch( incore_task_data& task, const basis_type& basis, const Molecule& mol, const value_type* P,
                           int64_t ldp, value_type* VXC, int64_t ldvxc, value_type* EXC,
                           value_type* N_EL, incore_integrator_type& incore_integrator);
public:

  template <typename... Args>
  ShellBatchedPGASDistributedXCIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ShellBatchedPGASDistributedXCIntegrator() noexcept = default;

};

}
}
