#pragma once


#include <gauxc/xc_integrator/integrator_defaults.hpp>

namespace GauXC {

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl ) :
  pimpl_( std::move( pimpl ) ) { }

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( MPI_Comm comm, const functional_type& func, 
  const basisset_type& basis, std::shared_ptr<LoadBalancer> lb) :
  XCIntegrator( detail::make_default_host_integrator<MatrixType>(
    comm,
    std::make_shared<functional_type>(func),
    std::make_shared<basisset_type>(basis),
    lb
  )) { }


template <typename MatrixType>
XCIntegrator<MatrixType>::~XCIntegrator() noexcept = default;



template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_vxc_type
  XCIntegrator<MatrixType>::eval_exc_vxc( const MatrixType& P ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");

  return pimpl_->eval_exc_vxc(P);
};
}
