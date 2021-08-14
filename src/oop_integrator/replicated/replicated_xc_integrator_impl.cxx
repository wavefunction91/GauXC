#include <gauxc/oop_xc_integrator/replicated/replicated_xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type > func,
                              std::shared_ptr< LoadBalancer >    lb, 
                              std::unique_ptr< LocalWorkDriver >&& lwd ) :
    func_(func), load_balancer_(lb), local_work_driver_(std::move(lwd)) { }

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ~ReplicatedXCIntegratorImpl() noexcept = default;


template <typename ValueType>
void ReplicatedXCIntegratorImpl<ValueType>::
  eval_exc_vxc( int64_t m, int64_t n, const value_type* P,
                int64_t ldp, value_type* VXC, int64_t ldvxc,
                value_type* EXC ) {

    eval_exc_vxc_(m,n,P,ldp,VXC,ldvxc,EXC);

}

template class ReplicatedXCIntegratorImpl<double>;

}
}
