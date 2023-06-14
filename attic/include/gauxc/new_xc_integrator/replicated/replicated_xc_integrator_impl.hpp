#pragma once

#include <gauxc/new_xc_integrator/replicated_xc_integrator.hpp>
#include <gauxc/types.hpp>
#include <gauxc/basisset.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
class ReplicatedXCIntegratorImpl {

public:

  using value_type = ValueType;
  using basis_type = BasisSet< value_type >;

protected:

#ifdef GAUXC_ENABLE_MPI
  MPI_Comm comm_;
#endif

  std::shared_ptr< functional_type > func_;
  std::shared_ptr< basis_type >      basis_;

  std::shared_ptr< LoadBalancer >    load_balancer_;

  util::Timer timer_;


  virtual void eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                              int64_t ldp, value_type* VXC, int64_t ldvxc,
                              value_type* EXC ) = 0;
public:

#ifdef GAUXC_ENABLE_MPI

  ReplicatedXCIntegratorImpl( MPI_Comm comm,
                              std::shared_ptr< functional_type > func,
                              std::shared_ptr< basis_type >      basis,
                              std::shared_ptr< LoadBalancer >    lb );

#else

  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type > func,
                              std::shared_ptr< basis_type >      basis,
                              std::shared_ptr< LoadBalancer >    lb );

#endif

  virtual ~ReplicatedXCIntegratorImpl() noexcept;

  void eval_exc_vxc( int64_t m, int64_t n, const value_type* P,
                     int64_t ldp, value_type* VXC, int64_t ldvxc,
                     value_type* EXC ); 

  inline const util::Timer& get_timings() const { return timer_; }

};


extern template class ReplicatedXCIntegratorImpl<double>;

}
}
