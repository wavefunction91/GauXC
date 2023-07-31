/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once
#include <gauxc/xc_integrator/replicated/replicated_xc_device_integrator.hpp>
#include "device/xc_device_data.hpp"
#include "incore_replicated_xc_device_integrator.hpp"

namespace GauXC {
namespace detail {

template <typename ValueType>
class ShellBatchedReplicatedXCDeviceIntegrator : 
  public ReplicatedXCDeviceIntegrator<ValueType> {

  using base_type  = ReplicatedXCDeviceIntegrator<ValueType>;

public:

  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

  using host_task_container = std::vector<XCTask>;
  using host_task_iterator  = typename host_task_container::iterator;

protected:

  using incore_integrator_type =
    IncoreReplicatedXCDeviceIntegrator<ValueType>;

  // Struct to manage data associated with task subset to execute on the device
  struct incore_device_task {
    host_task_iterator   task_begin;
    host_task_iterator   task_end;
    std::vector<int32_t> shell_list;
  };

  void integrate_den_( int64_t m, int64_t n, const value_type* P,
                       int64_t ldp, value_type* integrate_den ) override;

  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* VXC, int64_t ldvxc,
                      value_type* EXC ) override;

  void eval_exc_vxc_UKS_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp,
                      int64_t m1, int64_t n1, const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXC, int64_t ldvxc,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC ) override;

  void eval_exc_vxc_GKS_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp,
                      int64_t m1, int64_t n1, const value_type* Pz,
                      int64_t ldpz,
                      int64_t m2, int64_t n2, const value_type* Px,
                      int64_t ldpx,
                      int64_t m3, int64_t n3, const value_type* Py,
                      int64_t ldpy,
                      value_type* VXC, int64_t ldvxc,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* VXCx, int64_t ldvxcx,
                      value_type* VXCy, int64_t ldvxcy,
                      value_type* EXC ) override;

  void eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                       int64_t ldp, value_type* EXC_GRAD ) override;

  void eval_exx_( int64_t m, int64_t n, const value_type* P,
                  int64_t ldp, value_type* K, int64_t ldk,
                  const IntegratorSettingsEXX& settings ) override;

  void exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                            value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            incore_integrator_type& incore_integrator,
                            XCDeviceData& device_data );

  void eval_exc_grad_local_work_( int64_t m, int64_t n, const value_type* P,
                                  int64_t ldp, value_type* EXC_GRAD, 
                                  host_task_iterator task_begin, host_task_iterator task_end,
                                  incore_integrator_type& incore_integrator,
                                  XCDeviceData& device_data );

  
  incore_device_task generate_incore_device_task( const uint32_t     nbf_threshold,
                                                  const basis_type&  basis,
                                                  host_task_iterator task_begin,
                                                  host_task_iterator task_end );

  void execute_task_batch( incore_device_task& task, const basis_type& basis, const Molecule& mol, const value_type* P,
                           int64_t ldp, value_type* VXC, int64_t ldvxc, value_type* EXC,
                           value_type* N_EL, incore_integrator_type& incore_integrator,
                           XCDeviceData& device_data );
public:

  template <typename... Args>
  ShellBatchedReplicatedXCDeviceIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ShellBatchedReplicatedXCDeviceIntegrator() noexcept;

};

extern template class ShellBatchedReplicatedXCDeviceIntegrator<double>;

}
}
