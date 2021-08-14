#pragma once
#include <gauxc/oop_xc_integrator/replicated/replicated_xc_host_integrator.hpp>
#include "xc_host_data.hpp"
#include "integrator_util/xc_integrator_state.hpp"

namespace GauXC {
namespace detail {

template <typename ValueType>
class ReferenceReplicatedXCHostIntegrator : 
  public ReplicatedXCHostIntegrator<ValueType> {

  using base_type  = ReplicatedXCHostIntegrator<ValueType>;

public:

  using value_type = typename base_type::value_type;
  using basis_type = typename base_type::basis_type;

protected:

  XCIntegratorState state_;

  void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                      int64_t ldp, value_type* VXC, int64_t ldvxc,
                      value_type* EXC ) override;

  void exc_vxc_local_work_( const value_type* P, int64_t ldp, value_type* VXC,
                            int64_t ldvxc, value_type* EXC, value_type *N_EL,
                            XCHostData<value_type>& host_data );
public:

  template <typename... Args>
  ReferenceReplicatedXCHostIntegrator( Args&&... args ) :
    base_type( std::forward<Args>(args)... ) { }

  virtual ~ReferenceReplicatedXCHostIntegrator() noexcept;

};

extern template class ReferenceReplicatedXCHostIntegrator<double>;

}
}
