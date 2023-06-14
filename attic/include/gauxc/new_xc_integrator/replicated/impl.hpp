#pragma once

#include <gauxc/new_xc_integrator/replicated/replicated_xc_integrator_impl.hpp>

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
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  matrix_type VXC( P.rows(), P.cols() );
  value_type  EXC;

  if( not pimpl_ ) throw std::runtime_error( "Not Initialized" );
  pimpl_->eval_exc_vxc( P.rows(), P.cols(), P.data(), P.rows(),
                        VXC.data(), VXC.rows(), &EXC );

  return std::make_tuple( EXC, VXC );

}

}
}
