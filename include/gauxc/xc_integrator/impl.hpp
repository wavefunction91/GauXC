/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
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
typename XCIntegrator<MatrixType>::value_type
  XCIntegrator<MatrixType>::integrate_den( const MatrixType& P ) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->integrate_den(P);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_vxc_type_rks
  XCIntegrator<MatrixType>::eval_exc_vxc( const MatrixType& P ) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->eval_exc_vxc(P);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_vxc_type_uks
  XCIntegrator<MatrixType>::eval_exc_vxc( const MatrixType& Pscalar, const MatrixType& Pz ) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->eval_exc_vxc(Pscalar, Pz);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::exc_grad_type
  XCIntegrator<MatrixType>::eval_exc_grad( const MatrixType& P ) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->eval_exc_grad(P);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::exx_type
  XCIntegrator<MatrixType>::eval_exx( const MatrixType&     P,
                                      const IntegratorSettingsEXX& settings ) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->eval_exx(P,settings);
};

template <typename MatrixType>
typename XCIntegrator<MatrixType>::atomic_overlap_type
  XCIntegrator<MatrixType>::eval_atomic_overlap(int64_t iAtom) {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->eval_atomic_overlap(iAtom);
};

template <typename MatrixType>
const util::Timer& XCIntegrator<MatrixType>::get_timings() const {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_timings();
}

template <typename MatrixType>
const LoadBalancer& XCIntegrator<MatrixType>::load_balancer() const {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->load_balancer();
}

template <typename MatrixType>
LoadBalancer& XCIntegrator<MatrixType>::load_balancer() {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->load_balancer();
}
 
}
