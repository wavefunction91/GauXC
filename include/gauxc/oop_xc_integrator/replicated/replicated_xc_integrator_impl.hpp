#pragma once

#include <gauxc/oop_xc_integrator/replicated_xc_integrator.hpp>
#include <gauxc/oop_xc_integrator/local_work_driver.hpp>
#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>

namespace GauXC  {
namespace detail {


/// Base class for ReplicatedXCIntegrator implementations
template <typename ValueType>
class ReplicatedXCIntegratorImpl {

public:

  using value_type = ValueType;
  using basis_type = BasisSet< value_type >;

protected:

  std::shared_ptr< functional_type > func_;               ///< XC functional
  std::shared_ptr< LoadBalancer >    load_balancer_;      ///< Load Balancer
  std::unique_ptr< LocalWorkDriver > local_work_driver_;  ///< Local Work Driver

  util::Timer timer_;


  virtual void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                              int64_t ldp, value_type* VXC, int64_t ldvxc,
                              value_type* EXC ) = 0;
public:

  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type > func,
                              std::shared_ptr< LoadBalancer >    lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd );

  virtual ~ReplicatedXCIntegratorImpl() noexcept;

  void eval_exc_vxc( int64_t m, int64_t n, const value_type* P,
                     int64_t ldp, value_type* VXC, int64_t ldvxc,
                     value_type* EXC ); 

  inline const util::Timer& get_timings() const { return timer_; }

};


extern template class ReplicatedXCIntegratorImpl<double>;

}
}
