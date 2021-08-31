#pragma once

#include <gauxc/xc_integrator/xc_integrator_impl.hpp>

// Implementations of XCIntegrator public API

namespace GauXC {

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator( std::unique_ptr<pimpl_type>&& pimpl ) :
  pimpl_( std::move( pimpl ) ) { }


template <typename MatrixType>
XCIntegrator<MatrixType>::~XCIntegrator() noexcept = default;

template <typename MatrixType>
XCIntegrator<MatrixType>::XCIntegrator(XCIntegrator&&) noexcept = default;


template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_vxc_type
  XCIntegrator<MatrixType>::eval_exc_vxc( const MatrixType& P ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->eval_exc_vxc(P);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_grad_type
  XCIntegrator<MatrixType>::eval_exc_grad( const MatrixType& P ) {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");
  return pimpl_->eval_exc_grad(P);
};

template <typename MatrixType>
const util::Timer& XCIntegrator<MatrixType>::get_timings() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");

  return pimpl_->get_timings();
}

template <typename MatrixType>
const LoadBalancer& XCIntegrator<MatrixType>::load_balancer() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");

  return pimpl_->load_balancer();
}
 
}
