/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include <gauxc/xc_integrator/pgas_dist/pgas_dist_xc_integrator_impl.hpp>
#include <gauxc/exceptions.hpp>

// Implementations of PGASDistributedXCIntegrator public API

namespace GauXC  {
namespace detail {


template <typename MatrixType>
PGASDistributedXCIntegrator<MatrixType>::
  PGASDistributedXCIntegrator( std::unique_ptr<pimpl_type>&& pimpl ) : 
    pimpl_(std::move(pimpl)){ }

template <typename MatrixType>
PGASDistributedXCIntegrator<MatrixType>::PGASDistributedXCIntegrator(): 
  PGASDistributedXCIntegrator(nullptr){ }

template <typename MatrixType>
PGASDistributedXCIntegrator<MatrixType>::~PGASDistributedXCIntegrator() noexcept = default; 
template <typename MatrixType>
PGASDistributedXCIntegrator<MatrixType>::
  PGASDistributedXCIntegrator(PGASDistributedXCIntegrator&&) noexcept = default; 

template <typename MatrixType>
const util::Timer& PGASDistributedXCIntegrator<MatrixType>::get_timings_() const {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_timings();
}

template <typename MatrixType>
const LoadBalancer& PGASDistributedXCIntegrator<MatrixType>::get_load_balancer_() const {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_load_balancer();
}
template <typename MatrixType>
LoadBalancer& PGASDistributedXCIntegrator<MatrixType>::get_load_balancer_() {
  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  return pimpl_->get_load_balancer();
}


template <typename MatrixType>
typename PGASDistributedXCIntegrator<MatrixType>::value_type 
  PGASDistributedXCIntegrator<MatrixType>::integrate_den_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  value_type N_EL;
  
  //pimpl_->integrate_den( P.rows(), P.cols(), P.data(), P.rows(), &N_EL );

  return N_EL;
}

template <typename MatrixType>
typename PGASDistributedXCIntegrator<MatrixType>::exc_vxc_type_rks 
  PGASDistributedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  //matrix_type VXC( P.rows(), P.cols() );
  matrix_type VXC;
  value_type  EXC;
  if constexpr (std::is_convertible_v<MatrixType, Darray>)
    pimpl_->eval_exc_vxc( P, VXC, &EXC );
  else
    GAUXC_GENERIC_EXCEPTION("Passed density not compatible with internal PGAS matrix");

  return std::make_tuple( EXC, VXC );

}

template <typename MatrixType>
typename PGASDistributedXCIntegrator<MatrixType>::exc_vxc_type_uks
  PGASDistributedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& Pscalar, const MatrixType& Pz ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  //matrix_type VXCscalar( Pscalar.rows(), Pscalar.cols() );
  //matrix_type VXCz( Pz.rows(), Pz.cols() );
  matrix_type VXCscalar, VXCz;
  value_type  EXC;

  //pimpl_->eval_exc_vxc( Pscalar.rows(), Pscalar.cols(), Pscalar.data(), Pscalar.rows(),
  //                      Pz.data(), Pz.rows(),
  //                      VXCscalar.data(), VXCscalar.rows(),
  //                      VXCz.data(), VXCz.rows(), &EXC );

  return std::make_tuple( EXC, VXCscalar, VXCz );

}

template <typename MatrixType>
typename PGASDistributedXCIntegrator<MatrixType>::exc_grad_type 
  PGASDistributedXCIntegrator<MatrixType>::eval_exc_grad_( const MatrixType& P ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();

  std::vector<value_type> EXC_GRAD( 3*pimpl_->load_balancer().molecule().natoms() );
  //pimpl_->eval_exc_grad( P.rows(), P.cols(), P.data(), P.rows(),
  //                       EXC_GRAD.data() );

  return EXC_GRAD;

}

template <typename MatrixType>
typename PGASDistributedXCIntegrator<MatrixType>::exx_type 
  PGASDistributedXCIntegrator<MatrixType>::eval_exx_( const MatrixType& P, const IntegratorSettingsEXX& settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  
  matrix_type K;

  //pimpl_->eval_exx( P.rows(), P.cols(), P.data(), P.rows(),
  //                  K.data(), K.rows(), settings );

  return K;

}

}
}
