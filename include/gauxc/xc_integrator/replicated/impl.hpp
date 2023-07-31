/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/replicated/replicated_xc_integrator_impl.hpp>
#include <gauxc/exceptions.hpp>

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
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_timings();
}

template <typename MatrixType>
const LoadBalancer& ReplicatedXCIntegrator<MatrixType>::get_load_balancer_() const {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_load_balancer();
}
template <typename MatrixType>
LoadBalancer& ReplicatedXCIntegrator<MatrixType>::get_load_balancer_() {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_load_balancer();
}


template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::value_type 
  ReplicatedXCIntegrator<MatrixType>::integrate_den_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  value_type N_EL;
  
  pimpl_->integrate_den( P.rows(), P.cols(), P.data(), P.rows(), &N_EL );

  return N_EL;
}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXC( P.rows(), P.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc( P.rows(), P.cols(), P.data(), P.rows(),
                        VXC.data(), VXC.rows(), &EXC );

  return std::make_tuple( EXC, VXC );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_UKS
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_UKS_( const MatrixType& P, const MatrixType& Pz ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXC( P.rows(), P.cols() );
  matrix_type VXCz( Pz.rows(), Pz.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc_UKS( P.rows(), P.cols(), P.data(), P.rows(),
                        Pz.rows(), Pz.cols(), Pz.data(), Pz.rows(),
                        VXC.data(), VXC.rows(),
                        VXCz.data(), VXCz.rows(), &EXC );

  return std::make_tuple( EXC, VXC, VXCz );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_GKS
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_GKS_( const MatrixType& P,const MatrixType& Pz, const MatrixType& Px,const MatrixType& Py ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXC( P.rows(), P.cols() );
  matrix_type VXCz( Pz.rows(), Pz.cols() );
  matrix_type VXCx( Px.rows(), Px.cols() );
  matrix_type VXCy( Py.rows(), Py.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc_GKS( P.rows(), P.cols(), P.data(), P.rows(),
                        Pz.rows(), Pz.cols(), Pz.data(), Pz.rows(),
                        Px.rows(), Px.cols(), Px.data(), Px.rows(),
                        Py.rows(), Py.cols(), Py.data(), Py.rows(),
                        VXC.data(), VXC.rows(),
                        VXCz.data(), VXCz.rows(),
                        VXCx.data(), VXCx.rows(),
                        VXCy.data(), VXCy.rows(), &EXC );

  return std::make_tuple( EXC, VXC, VXCz, VXCx, VXCy );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_grad_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_grad_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();

  std::vector<value_type> EXC_GRAD( 3*pimpl_->load_balancer().molecule().natoms() );
  pimpl_->eval_exc_grad( P.rows(), P.cols(), P.data(), P.rows(),
                         EXC_GRAD.data() );

  return EXC_GRAD;

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exx_type 
  ReplicatedXCIntegrator<MatrixType>::eval_exx_( const MatrixType& P, const IntegratorSettingsEXX& settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  
  matrix_type K( P.rows(), P.cols() );

  pimpl_->eval_exx( P.rows(), P.cols(), P.data(), P.rows(),
                    K.data(), K.rows(), settings );

  return K;

}

}
}
