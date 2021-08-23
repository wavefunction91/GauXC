#pragma once
#include <gauxc/oop_xc_integrator/replicated/replicated_xc_device_integrator.hpp>
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

  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* VXC, int64_t ldvxc,
                      value_type* EXC ) override;

  void exc_vxc_local_work_( const basis_type& basis, const value_type* P, int64_t ldp, 
                            value_type* VXC, int64_t ldvxc, value_type* EXC, value_type *N_EL,
                            host_task_iterator task_begin, host_task_iterator task_end,
                            incore_integrator_type& incore_integrator,
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
