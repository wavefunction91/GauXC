/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
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
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_rks 
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& P, const IntegratorSettingsXC& ks_settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXC( P.rows(), P.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc( P.rows(), P.cols(), P.data(), P.rows(),
                        VXC.data(), VXC.rows(), &EXC, ks_settings );

  return std::make_tuple( EXC, VXC );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_uks
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& Ps, const MatrixType& Pz, const IntegratorSettingsXC& ks_settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXCs( Ps.rows(), Ps.cols() );
  matrix_type VXCz( Pz.rows(), Pz.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc( Ps.rows(), Ps.cols(), Ps.data(), Ps.rows(),
                        Pz.data(), Pz.rows(),
                        VXCs.data(), VXCs.rows(),
                        VXCz.data(), VXCz.rows(), &EXC, ks_settings );

  return std::make_tuple( EXC, VXCs, VXCz );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_gks
  ReplicatedXCIntegrator<MatrixType>::eval_exc_vxc_( const MatrixType& Ps, const MatrixType& Pz, const MatrixType& Py, const MatrixType& Px,
                                                     const IntegratorSettingsXC& ks_settings) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type VXCs( Ps.rows(), Ps.cols() );
  matrix_type VXCz( Pz.rows(), Pz.cols() );
  matrix_type VXCy( Py.rows(), Py.cols() );
  matrix_type VXCx( Px.rows(), Px.cols() );
  value_type  EXC;

  pimpl_->eval_exc_vxc( Ps.rows(), Ps.cols(), Ps.data(), Ps.rows(),
                        Pz.data(), Pz.rows(),
                        Py.data(), Py.rows(),
                        Px.data(), Px.rows(),
                        VXCs.data(), VXCs.rows(),
                        VXCz.data(), VXCz.rows(),
                        VXCy.data(), VXCy.rows(),
                        VXCx.data(), VXCx.rows(), &EXC, ks_settings );

  return std::make_tuple( EXC, VXCs, VXCz, VXCy, VXCx);

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_neo_rks
  ReplicatedXCIntegrator<MatrixType>::neo_eval_exc_vxc_( const MatrixType& elec_Ps, const MatrixType& prot_Ps, const MatrixType& prot_Pz,
                                                         const IntegratorSettingsXC& ks_settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type elec_VXCs( elec_Ps.rows(), elec_Ps.cols() );
  matrix_type prot_VXCs( prot_Ps.rows(), prot_Ps.cols() );
  matrix_type prot_VXCz( prot_Pz.rows(), prot_Pz.cols() );
  value_type  EXC;

  pimpl_->neo_eval_exc_vxc( elec_Ps.rows(), elec_Ps.cols(), prot_Ps.rows(), prot_Ps.cols(),
                            elec_Ps.data(), elec_Ps.rows(),
                            prot_Ps.data(), prot_Ps.rows(),
                            prot_Pz.data(), prot_Pz.rows(),
                            elec_VXCs.data(), elec_VXCs.rows(),
                            prot_VXCs.data(), prot_VXCs.rows(),
                            prot_VXCz.data(), prot_VXCz.rows(),
                            &EXC);

  return std::make_tuple( EXC, elec_VXCs, prot_VXCs, prot_VXCz );

}

template <typename MatrixType>
typename ReplicatedXCIntegrator<MatrixType>::exc_vxc_type_neo_uks
  ReplicatedXCIntegrator<MatrixType>::neo_eval_exc_vxc_( const MatrixType& elec_Ps, const MatrixType& elec_Pz, const MatrixType& prot_Ps, const MatrixType& prot_Pz,
                                                         const IntegratorSettingsXC& ks_settings ) {

  if( not pimpl_ ) GAUXC_PIMPL_NOT_INITIALIZED();
  matrix_type elec_VXCs( elec_Ps.rows(), elec_Ps.cols() );
  matrix_type elec_VXCz( elec_Pz.rows(), elec_Pz.cols() );
  matrix_type prot_VXCs( prot_Ps.rows(), prot_Ps.cols() );
  matrix_type prot_VXCz( prot_Pz.rows(), prot_Pz.cols() );
  value_type  EXC;

  pimpl_->neo_eval_exc_vxc( elec_Ps.rows(), elec_Ps.cols(), prot_Ps.rows(), prot_Ps.cols(),
                            elec_Ps.data(), elec_Ps.rows(),
                            elec_Pz.data(), elec_Pz.rows(),
                            prot_Ps.data(), prot_Ps.rows(),
                            prot_Pz.data(), prot_Pz.rows(),
                            elec_VXCs.data(), elec_VXCs.rows(),
                            elec_VXCz.data(), elec_VXCz.rows(),
                            prot_VXCs.data(), prot_VXCs.rows(),
                            prot_VXCz.data(), prot_VXCz.rows(),
                            &EXC);

  return std::make_tuple( EXC, elec_VXCs, elec_VXCz, prot_VXCs, prot_VXCz );

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
