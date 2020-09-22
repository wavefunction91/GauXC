#pragma once

#include <gauxc/xc_integrator/integrator_factory.hpp>

namespace GauXC {

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl ) :
  pimpl_( std::move( pimpl ) ) { }


#ifdef GAUXC_ENABLE_MPI

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( ExecutionSpace ex, MPI_Comm comm, 
                                        const functional_type& func, 
                                        const basisset_type& basis, 
                                        std::shared_ptr<LoadBalancer> lb) :
  XCIntegrator( detail::integrator_factory<MatrixType>(
    ex, comm, func, basis, lb
  )) { }

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( MPI_Comm comm, const functional_type& func, 
                                        const basisset_type& basis, 
                                        std::shared_ptr<LoadBalancer> lb) :
  XCIntegrator( ExecutionSpace::Host, comm, func, basis, lb ) { };

#else

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( ExecutionSpace ex, 
                                        const functional_type& func, 
                                        const basisset_type& basis, 
                                        std::shared_ptr<LoadBalancer> lb) :
  XCIntegrator( detail::integrator_factory<MatrixType>(
    ex, func, basis, lb
  )) { }

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( const functional_type& func, 
                                        const basisset_type& basis, 
                                        std::shared_ptr<LoadBalancer> lb) :
  XCIntegrator( ExecutionSpace::Host, func, basis, lb ) { };

#endif

template <typename MatrixType>
XCIntegrator<MatrixType>::~XCIntegrator() noexcept = default;



template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_vxc_type
  XCIntegrator<MatrixType>::eval_exc_vxc( const MatrixType& P ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");

  return pimpl_->eval_exc_vxc(P);
};
}
