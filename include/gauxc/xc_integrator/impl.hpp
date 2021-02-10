#pragma once

#include <gauxc/xc_integrator/integrator_factory.hpp>

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
const util::Timer& XCIntegrator<MatrixType>::get_timings() const {
  if( not pimpl_ ) throw std::runtime_error("Not Initialized");

  return pimpl_->get_timings();
}
}
