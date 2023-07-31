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

namespace GauXC {
namespace detail {

template <typename ValueType>
class IncoreReplicatedXCDeviceIntegrator : 
  public ReplicatedXCDeviceIntegrator<ValueType> {

  using base_type  = ReplicatedXCDeviceIntegrator<ValueType>;

public:

  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

  using host_task_container = std::vector<XCTask>;
  using host_task_iterator  = typename host_task_container::iterator;

protected:

  void integrate_den_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* N_EL ) override;

  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* VXC, int64_t ldvxc,
                      value_type* EXC ) override;

  void eval_exc_vxc_UKS_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp,
                      const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXC, int64_t ldvxc,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC ) override;

  void eval_exc_vxc_GKS_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp,
                      const value_type* Pz,
                      int64_t ldpz,
                      const value_type* Px,
                      int64_t ldpx,
                      const value_type* Py,
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


  void integrate_den_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                            value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );

  void exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );

  void exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                            value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );

  void exc_vxc_local_work_UKS_( const basis_type& basis, const value_type* P, int64_t ldp,
                                const value_type* Pz, int64_t ldpz,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );

  void exc_vxc_local_work_UKS_( const basis_type& basis, const value_type* P, int64_t ldp,
                            const value_type* Pz, int64_t ldpz,
                            value_type* VXC, int64_t ldvxc,
                            value_type* VXCz, int64_t ldvxcz, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );


  void exc_vxc_local_work_GKS_( const basis_type& basis, const value_type* P, int64_t ldp,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Px, int64_t ldpx,
                            const value_type* Py, int64_t ldpy,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );

  void exc_vxc_local_work_GKS_( const basis_type& basis, const value_type* P, int64_t ldp,
                            const value_type* Pz, int64_t ldpz,
                            const value_type* Px, int64_t ldpx,
                            const value_type* Py, int64_t ldpy,
                            value_type* VXC, int64_t ldvxc,
                            value_type* VXCz, int64_t ldvxcz,
                            value_type* VXCx, int64_t ldvxcx,
                            value_type* VXCy, int64_t ldvxcy, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            XCDeviceData& device_data );





  void eval_exc_grad_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                                  host_task_iterator task_begin, host_task_iterator task_end,
                                  XCDeviceData& device_data );

  void eval_exc_grad_local_work_( const basis_type& basis, const value_type* P,
                                  int64_t ldp, value_type* EXC_GRAD, 
                                  host_task_iterator task_begin, host_task_iterator task_end,
                                  XCDeviceData& device_data );



  void exx_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                        host_task_iterator task_begin, host_task_iterator task_end,
                        XCDeviceData& device_data, 
                        const IntegratorSettingsEXX& settings);

  void exx_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                        value_type* K, int64_t ldk,
                        host_task_iterator task_begin, host_task_iterator task_end,
                        XCDeviceData& device_data, 
                        const IntegratorSettingsEXX& settings);

  void exx_ek_screening_local_work_( const basis_type& basis, 
                        const value_type* P, int64_t ldp, 
                        XCDeviceData& device_data, 
                        const IntegratorSettingsEXX& settings);

public:

  template <typename... Args>
  IncoreReplicatedXCDeviceIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~IncoreReplicatedXCDeviceIntegrator() noexcept;


  template <typename... Args>
  void exc_vxc_local_work(Args&&... args) {
    exc_vxc_local_work_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void exc_vxc_local_work_UKS(Args&&... args) {
    exc_vxc_local_work_UKS_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void exc_vxc_local_work_GKS(Args&&... args) {
    exc_vxc_local_work_GKS_( std::forward<Args>(args)... );
  }

};

extern template class IncoreReplicatedXCDeviceIntegrator<double>;

}
}
