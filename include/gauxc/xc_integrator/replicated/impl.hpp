#pragma once

#include <gauxc/xc_integrator/replicated/replicated_xc_integrator_impl.hpp>

// Implementations of ReplicatedXCIntegrator public API

namespace GauXC  {
namespace detail {


template <typename MatrixType>
ReplicatedXCIntegrator<MatrixType>::
  ReplicatedXCIntegrator( std::unique_ptr<pimpl_type>&& pimpl ) : 
    pimpl_(std::move(pimpl)){ }

template <typename MatrixType>
ReplicatedXCIntegrator<MatrixType>::ReplicatedXCIntegrator(): 
  ReplicatedXCIntegrator(nullptr){ }

template <typename MatrixType>
ReplicatedXCIntegrator<MatrixType>::~ReplicatedXCIntegrator() noexcept = default; 
template <typename MatrixType>
ReplicatedXCIntegrator<MatrixType>::
  ReplicatedXCIntegrator(ReplicatedXCIntegrator&&) noexcept = default; 

template <typename MatrixType>
const util::Timer& ReplicatedXCIntegrator<MatrixType>::get_timings_() const {
  if( not pimpl_ ) throw std::runtime_error( "Not Initialized" );
  return pimpl_->get_timings();
}

template <typename MatrixType>
const LoadBalancer& ReplicatedXCIntegrator<MatrixType>::get_load_balancer_() const {
  if( not pimpl_ ) throw std::runtime_error( "Not Initialized" );
  return pimpl_->get_load_balancer();
}


template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  matrix_type VXC( P.rows(), P.cols() );
  value_type  EXC;

  if( not pimpl_ ) throw std::runtime_error( "Not Initialized" );
  pimpl_->eval_exc_vxc( P.rows(), P.cols(), P.data(), P.rows(),
                        VXC.data(), VXC.rows(), &EXC );

  return std::make_tuple( EXC, VXC );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_grad_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_grad_( const MatrixType& P ) {

  if( not pimpl_ ) throw std::runtime_error( "Not Initialized" );

  std::vector<value_type> EXC_GRAD( 3*pimpl_->load_balancer().molecule().natoms() );
  pimpl_->eval_exc_grad( P.rows(), P.cols(), P.data(), P.rows(),
                         EXC_GRAD.data() );

  return EXC_GRAD;

}

}
}
