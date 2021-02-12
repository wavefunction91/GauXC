#include <gauxc/new_xc_integrator/replicated/replicated_xc_integrator_impl.hpp>

namespace GauXC  {
namespace detail {

#ifdef GAUXC_ENABLE_MPI

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ReplicatedXCIntegratorImpl( MPI_Comm comm,
                              std::shared_ptr< functional_type > func,
                              std::shared_ptr< basis_type >      basis,
                              std::shared_ptr< LoadBalancer >    lb ) :
    comm_(comm), func_(func), basis_(basis), load_balancer_(lb) { }

#else

template <typename ValueType>
ReplicatedXCIntegratorImpl<ValueType>::
  ReplicatedXCIntegratorImpl( std::shared_ptr< functional_type > func,
                              std::shared_ptr< basis_type >      basis,
                              std::shared_ptr< LoadBalancer >    lb ) :
    func_(func), basis_(basis), load_balancer_(lb) { }

#endif

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
