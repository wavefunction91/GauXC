/**
 * GauXC Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/local_host_work_driver.hpp"
#include "host/blas.hpp"
#include <stdexcept>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC, const IntegratorSettingsXC& ks_settings ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxc < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");


  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    //exc_vxc_local_work_( P, ldp, VXC, ldvxc, EXC, &N_EL );
    exc_vxc_local_work_( P, ldp, nullptr, 0, nullptr, 0, nullptr, 0,
                         VXC, ldvxc, nullptr, 0, nullptr, 0, nullptr, 0, EXC, &N_EL, ks_settings );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXC, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );

  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldps < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPSCALAR");
  if( ldpz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");
  if( ldvxcs < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCSCALAR");
  if( ldvxcz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCZ");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( Ps, ldps, Pz, ldpz, nullptr, 0,nullptr, 0,
                         VXCs, ldvxcs, VXCz, ldvxcz, nullptr, 0, nullptr, 0, EXC, &N_EL, ks_settings );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );

  });


}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* Ps,
                      int64_t ldps,
                      const value_type* Pz,
                      int64_t ldpz,
                      const value_type* Py,
                      int64_t ldpy,
                      const value_type* Px,
                      int64_t ldpx,
                      value_type* VXCs, int64_t ldvxcs,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* VXCy, int64_t ldvxcy,
                      value_type* VXCx, int64_t ldvxcx,
                      value_type* EXC, const IntegratorSettingsXC& ks_settings ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldps < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPSCALAR");
  if( ldpz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");
  if( ldpy < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPX");
  if( ldpx < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPY");
  if( ldvxcs < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCSCALAR");
  if( ldvxcz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCZ");
  if( ldvxcy < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCX");
  if( ldvxcx < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCY");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;
   
  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( Ps, ldps, Pz, ldpz, Py, ldpy, Px, ldpx, 
                         VXCs, ldvxcs, VXCz, ldvxcz,
                         VXCy, ldvxcy, VXCx, ldvxcx, EXC, &N_EL, ks_settings );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXCs, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXCy, nbf*nbf, ReductionOp::Sum ); 
    this->reduction_driver_->allreduce_inplace( VXCx, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );

  });


}           
  
template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  neo_eval_exc_vxc_( int64_t m1, int64_t n1, int64_t m2, int64_t n2, 
                     const value_type* P1s, int64_t ldp1s,
                     const value_type* P2s, int64_t ldp2s,
                     const value_type* P2z, int64_t ldp2z,
                     value_type* VXC1s, int64_t ldvxc1s,
                     value_type* VXC2s, int64_t ldvxc2s,
                     value_type* VXC2z, int64_t ldvxc2z,
                     value_type* EXC1,  value_type* EXC2, const IntegratorSettingsXC& ks_settings ){
  
  const auto& basis  = this->load_balancer_->basis();
  const auto& basis2 = this->load_balancer_->basis2();

  // Check that P / VXC are sane
  const int64_t nbf1 = basis.nbf();
  const int64_t nbf2 = basis2.nbf();

  if( m1 != n1 | m2 != n2)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m1 != nbf1 | m2 != nbf2)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp1s < nbf1 | ldp2s < nbf2 | ldp2z < nbf2 )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxc1s < nbf1 | ldvxc2s < nbf2 | ldvxc2z < nbf2 )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    neo_exc_vxc_local_work_( P1s, ldp1s,
                             nullptr, 0,
                             P2s, ldp2s,
                             P2z, ldp2z,
                             VXC1s, ldvxc1s,
                             nullptr, 0,
                             VXC2s, ldvxc2s,
                             VXC2z, ldvxc2z,
                             EXC1, EXC2, &N_EL, ks_settings  );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXC1s, nbf1*nbf1,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXC2s, nbf2*nbf1,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXC2z, nbf2*nbf2,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC1,  1    ,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC2,  1    ,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    ,  ReductionOp::Sum );

  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_neo_exc_vxc_( int64_t m1, int64_t n1, int64_t m2, int64_t n2,
                     const value_type* P1s, int64_t ldp1s,
                     const value_type* P1z, int64_t ldp1z,
                     const value_type* P2s, int64_t ldp2s,
                     const value_type* P2z, int64_t ldp2z,
                     value_type* VXC1s, int64_t ldvxc1s,
                     value_type* VXC1z, int64_t ldvxc1z,
                     value_type* VXC2s, int64_t ldvxc2s,
                     value_type* VXC2z, int64_t ldvxc2z,
                     value_type* EXC1, value_type* EXC2 ) {
  
  const auto& basis  = this->load_balancer_->basis();
  const auto& basis2 = this->load_balancer_->basis2();

  // Check that P / VXC are sane
  const int64_t nbf1 = basis.nbf();
  const int64_t nbf2 = basis2.nbf();

  if( m1 != n1 | m2 != n2)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m1 != nbf1 | m2 != nbf2)
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldp1s < nbf1 | ldp2s < nbf2 | ldp2z < nbf2 )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
  if( ldvxc1s < nbf1 | ldvxc2s < nbf2 | ldvxc2z < nbf2 )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXC");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    neo_exc_vxc_local_work_( P1s, ldp1s,
                             P1z, ldp1z,
                             P2s, ldp2s,
                             P2z, ldp2z,
                             VXC1s, ldvxc1s,
                             VXC1z, ldvxc1z,
                             VXC2s, ldvxc2s,
                             VXC2z, ldvxc2z,
                             EXC1, EXC2, &N_EL );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXC1s, nbf1*nbf1,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXC1z, nbf1*nbf1,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXC2s, nbf2*nbf1,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXC2z, nbf2*nbf2,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC1,  1    ,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC2,  1    ,  ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    ,  ReductionOp::Sum );

  });

}


template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_vxc_local_work_( const value_type* Ps, int64_t ldps,
                       const value_type* Pz, int64_t ldpz,
                       const value_type* Py, int64_t ldpy,
                       const value_type* Px, int64_t ldpx,
                       value_type* VXCs, int64_t ldvxcs,
                       value_type* VXCz, int64_t ldvxcz,
                       value_type* VXCy, int64_t ldvxcy,
                       value_type* VXCx, int64_t ldvxcx,
                       value_type* EXC, value_type *N_EL, const IntegratorSettingsXC& settings) {

  const bool is_gks = (Pz != nullptr) and (VXCz != nullptr) and (VXCy != nullptr) and (VXCx != nullptr);
  const bool is_uks = (Pz != nullptr) and (VXCz != nullptr) and (VXCy == nullptr) and (VXCx == nullptr);
  const bool is_rks = not is_uks and not is_gks;
  if (not is_rks and not is_uks and not is_gks) {
    GAUXC_GENERIC_EXCEPTION("MUST BE EITHER RKS, UKS, or GKS!");
  }


  // Misc KS settings
  IntegratorSettingsKS ks_settings;
  if( auto* tmp = dynamic_cast<const IntegratorSettingsKS*>(&settings) ) {
    ks_settings = *tmp;
  }

  const double gks_dtol = ks_settings.gks_dtol;

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();

  const bool needs_laplacian = func.needs_laplacian(); 
  
  if (func.is_mgga() and is_gks) {
    GAUXC_GENERIC_EXCEPTION("GKS NOT YET IMPLEMENTED WITH mGGA FUNCTIONALS!");
  }

  // Get basis map
  BasisSetMap basis_map(basis,mol);

  const int32_t nbf = basis.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };

  auto& tasks = this->load_balancer_->get_tasks();
  std::sort( tasks.begin(), tasks.end(), task_comparator );


  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified");
  }

  // Zero out integrands
  
  for( auto j = 0; j < nbf; ++j ) {
    for( auto i = 0; i < nbf; ++i ) {
      VXCs[i + j*ldvxcs] = 0.;
    }
  }
  if(not is_rks) {
    for( auto j = 0; j < nbf; ++j ) {
      for( auto i = 0; i < nbf; ++i ) {
        VXCz[i + j*ldvxcz] = 0.;
      }
    }
  }
  if(is_gks) {
    for( auto j = 0; j < nbf; ++j ) {
      for( auto i = 0; i < nbf; ++i ) {
        VXCy[i + j*ldvxcy] = 0.;
        VXCx[i + j*ldvxcx] = 0.;
      }
    }
  }
 
  double EXC_WORK = 0.0;
  double NEL_WORK = 0.0;
    
  // Loop over tasks
  const size_t ntasks = tasks.size();

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {
     
    //std::cout << iT << "/" << ntasks << std::endl;
    //printf("%lu / %lu\n", iT, ntasks);
    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate enough memory for batch
   
    const size_t spin_dim_scal = is_rks ? 1 : is_uks ? 2 : 4; // last case is_gks
    const size_t sds          = is_rks ? 1 : 2;
    const size_t gks_mod_KH = is_gks ? 6*npts : 0; // used to store H and H
    const size_t mgga_dim_scal = func.is_mgga() ? 4 : 1; // basis + d1basis

    // Things that every calc needs
    host_data.nbe_scr .resize(nbe  * nbe);
    host_data.zmat    .resize(npts * nbe * spin_dim_scal * mgga_dim_scal + gks_mod_KH); 
    host_data.eps     .resize(npts);
    host_data.vrho    .resize(npts * spin_dim_scal);

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts * spin_dim_scal);
    }
     
    // GGA data requirements
    const size_t gga_dim_scal = is_rks ? 1 : 3;
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( spin_dim_scal * 4 * npts );
      host_data.gamma      .resize( gga_dim_scal * npts );
      host_data.vgamma     .resize( gga_dim_scal * npts );
    }

    if( func.is_mgga() ){
      if ( needs_laplacian ) {
        host_data.basis_eval .resize( 11 * npts * nbe ); // basis + grad (3) + hess (6) + lapl 
        host_data.lapl       .resize( spin_dim_scal * npts );
        host_data.vlapl      .resize( spin_dim_scal * npts );
      } else {
        host_data.basis_eval .resize( 4 * npts * nbe ); // basis + grad (3)
      }

      host_data.den_scr    .resize( spin_dim_scal * 4 * npts );
      host_data.gamma      .resize( gga_dim_scal * npts );
      host_data.vgamma     .resize( gga_dim_scal * npts );
      host_data.tau        .resize( npts * spin_dim_scal );
      host_data.vtau       .resize( npts * spin_dim_scal );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    decltype(zmat) zmat_z = nullptr;
    decltype(zmat) zmat_x = nullptr;
    decltype(zmat) zmat_y = nullptr;
    if(!is_rks) {
      zmat_z = zmat + mgga_dim_scal * nbe * npts;
    }
    if(is_gks) {
      zmat_x = zmat_z + nbe * npts;
      zmat_y = zmat_x + nbe * npts;
    }
     
    auto* eps        = host_data.eps.data();
    auto* gamma      = host_data.gamma.data();
    auto* tau        = host_data.tau.data();
    auto* lapl       = host_data.lapl.data();
    auto* vrho       = host_data.vrho.data();
    auto* vgamma     = host_data.vgamma.data();
    auto* vtau       = host_data.vtau.data();
    auto* vlapl      = host_data.vlapl.data();


    value_type* dbasis_x_eval = nullptr;
    value_type* dbasis_y_eval = nullptr;
    value_type* dbasis_z_eval = nullptr;
    value_type* d2basis_xx_eval = nullptr;
    value_type* d2basis_xy_eval = nullptr;
    value_type* d2basis_xz_eval = nullptr;
    value_type* d2basis_yy_eval = nullptr;
    value_type* d2basis_yz_eval = nullptr;
    value_type* d2basis_zz_eval = nullptr;
    value_type* lbasis_eval = nullptr;
    value_type* dden_x_eval = nullptr;
    value_type* dden_y_eval = nullptr;
    value_type* dden_z_eval = nullptr;
    value_type* K = nullptr;
    value_type* H = nullptr;
    if (is_gks) { K = zmat + npts * nbe * 4; }
    value_type* mmat_x      = nullptr;
    value_type* mmat_y      = nullptr;
    value_type* mmat_z      = nullptr;
    value_type* mmat_x_z    = nullptr;
    value_type* mmat_y_z    = nullptr;
    value_type* mmat_z_z    = nullptr;

    if( func.is_gga() ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      dden_x_eval   = den_eval    + spin_dim_scal * npts;
      dden_y_eval   = dden_x_eval + spin_dim_scal * npts;
      dden_z_eval   = dden_y_eval + spin_dim_scal * npts;
      if (is_gks) { H = K + 3*npts;}
    }

    if ( func.is_mgga() ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      dden_x_eval   = den_eval    + spin_dim_scal * npts;
      dden_y_eval   = dden_x_eval + spin_dim_scal * npts;
      dden_z_eval   = dden_y_eval + spin_dim_scal * npts;
      mmat_x        = zmat + npts * nbe;
      mmat_y        = mmat_x + npts * nbe;
      mmat_z        = mmat_y + npts * nbe;
      if ( needs_laplacian ) {
        d2basis_xx_eval = dbasis_z_eval + npts * nbe;
        d2basis_xy_eval = d2basis_xx_eval + npts * nbe;
        d2basis_xz_eval = d2basis_xy_eval + npts * nbe;
        d2basis_yy_eval = d2basis_xz_eval + npts * nbe;
        d2basis_yz_eval = d2basis_yy_eval + npts * nbe;
        d2basis_zz_eval = d2basis_yz_eval + npts * nbe;
        lbasis_eval     = d2basis_zz_eval + npts * nbe;
      }
      if(is_uks) {
        mmat_x_z = zmat_z + npts * nbe;
        mmat_y_z = mmat_x_z + npts * nbe;
        mmat_z_z = mmat_y_z + npts * nbe;
      }
    }


    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation (+ Grad and Hessian)
    if( func.is_mgga() ) {
      if ( needs_laplacian ) {
        // TODO: Modify gau2grid to compute Laplacian instead of full hessian
        lwd->eval_collocation_hessian( npts, nshells, nbe, points, basis, shell_list,
          basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval,
          d2basis_xy_eval, d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval,
          d2basis_zz_eval);
        blas::lacpy( 'A', nbe, npts, d2basis_xx_eval, nbe, lbasis_eval, nbe );
        blas::axpy( nbe * npts, 1., d2basis_yy_eval, 1, lbasis_eval, 1);
        blas::axpy( nbe * npts, 1., d2basis_zz_eval, 1, lbasis_eval, 1);
      } else {
        lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
          basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
      }
    }
    // Evaluate Collocation (+ Grad)
    else if( func.is_gga() )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list,
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list,
        basis_eval );

     
    // Evaluate X matrix (fac * P * B) -> store in Z
    const auto xmat_fac = is_rks ? 2.0 : 1.0; // TODO Fix for spinor RKS input
    lwd->eval_xmat( mgga_dim_scal * npts, nbf, nbe, submat_map, xmat_fac, Ps, ldps, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    // X matrix for Pz
    if(not is_rks) {
      lwd->eval_xmat( mgga_dim_scal * npts, nbf, nbe, submat_map, 1.0, Pz, ldpz, basis_eval, nbe,
        zmat_z, nbe, nbe_scr);
    }
     
    if(is_gks) {
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, Py, ldpy, basis_eval, nbe,
        zmat_x, nbe, nbe_scr);
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, Px, ldpx, basis_eval, nbe,
        zmat_y, nbe, nbe_scr);
    }
     
    // Evaluate U and V variables
    if( func.is_mgga() ) {
      if (is_rks) {
        lwd->eval_uvvar_mgga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, lbasis_eval, zmat, nbe, mmat_x, mmat_y, mmat_z, 
          nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval, gamma, tau, lapl);
      } else if (is_uks) {
        lwd->eval_uvvar_mgga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, lbasis_eval, zmat, nbe, zmat_z, nbe, 
          mmat_x, mmat_y, mmat_z, nbe, mmat_x_z, mmat_y_z, mmat_z_z, nbe, 
          den_eval, dden_x_eval, dden_y_eval, dden_z_eval, gamma, tau, lapl);
      }
    } else if ( func.is_gga() ) {
      if(is_rks) {
        lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
          gamma );
      } else if(is_uks) {
        lwd->eval_uvvar_gga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, zmat_z, nbe, den_eval, dden_x_eval, 
          dden_y_eval, dden_z_eval, gamma );
      } else if(is_gks) {
        lwd->eval_uvvar_gga_gks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, zmat_z, nbe, zmat_x, nbe, zmat_y, nbe, den_eval, dden_x_eval,
          dden_y_eval, dden_z_eval, gamma, K, H, gks_dtol );
      }
       
     } else {
      if(is_rks) {
        lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, zmat, nbe, den_eval );
      } else if(is_uks) {
        lwd->eval_uvvar_lda_uks( npts, nbe, basis_eval, zmat, nbe, zmat_z, nbe,
          den_eval );
      } else if(is_gks) {
        lwd->eval_uvvar_lda_gks( npts, nbe, basis_eval, zmat, nbe, zmat_z, nbe,
          zmat_x, nbe, zmat_y, nbe, den_eval, K, gks_dtol );
      }
     }
    
    // Evaluate XC functional
    if( func.is_mgga() )
      func.eval_exc_vxc( npts, den_eval, gamma, lapl, tau, eps, vrho, vgamma, vlapl, vtau);
    else if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      vrho[sds*i] *= weights[i];
      if(not is_rks) vrho[sds*i+1] *= weights[i];
    }
     
    if( func.is_gga() ){
      for( int32_t i = 0; i < npts; ++i ) {
         vgamma[gga_dim_scal*i] *= weights[i];
         if(not is_rks) {
           vgamma[gga_dim_scal*i+1] *= weights[i];
           vgamma[gga_dim_scal*i+2] *= weights[i];
         }
      }
    }

    if( func.is_mgga() ){
      for( int32_t i = 0; i < npts; ++i) {
        vtau[spin_dim_scal*i]  *= weights[i];
        vgamma[gga_dim_scal*i] *= weights[i];
        if(not is_rks) {
          vgamma[gga_dim_scal*i+1] *= weights[i];
          vgamma[gga_dim_scal*i+2] *= weights[i];
          vtau[spin_dim_scal*i+1]  *= weights[i];
        }

        // TODO: Add checks for Lapacian-dependent functionals
        if( needs_laplacian ) {
          vlapl[spin_dim_scal*i] *= weights[i];
          if(not is_rks) {
            vlapl[spin_dim_scal*i+1] *= weights[i];
          }
        }
      }
    }



    // Evaluate Z matrix for VXC
    if( func.is_mgga() ) {
      if(is_rks) {
        lwd->eval_zmat_mgga_vxc_rks( npts, nbe, vrho, vgamma, vlapl, basis_eval, dbasis_x_eval,
                                     dbasis_y_eval, dbasis_z_eval, lbasis_eval,
                                     dden_x_eval, dden_y_eval, dden_z_eval, zmat, nbe);
        lwd->eval_mmat_mgga_vxc_rks( npts, nbe, vtau, vlapl, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval,
                                     mmat_x, mmat_y, mmat_z, nbe);
      } else if (is_uks) {
        lwd->eval_zmat_mgga_vxc_uks( npts, nbe, vrho, vgamma, vlapl, basis_eval, dbasis_x_eval,
                                     dbasis_y_eval, dbasis_z_eval, lbasis_eval,
                                     dden_x_eval, dden_y_eval, dden_z_eval, zmat, nbe, zmat_z, nbe);
        lwd->eval_mmat_mgga_vxc_uks( npts, nbe, vtau, vlapl, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval,
                                     mmat_x, mmat_y, mmat_z, nbe, mmat_x_z, mmat_y_z, mmat_z_z, nbe);
      }
    }
    else if( func.is_gga() ) {
      if(is_rks) {
        lwd->eval_zmat_gga_vxc_rks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe);
      } else if(is_uks) {
        lwd->eval_zmat_gga_vxc_uks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe, zmat_z, nbe);
      } else if(is_gks) {
        lwd->eval_zmat_gga_vxc_gks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe, zmat_z, nbe, zmat_x, nbe, zmat_y, nbe,
                                K, H);
      }
       
    } else {
      if(is_rks) {
        lwd->eval_zmat_lda_vxc_rks( npts, nbe, vrho, basis_eval, zmat, nbe );
      } else if(is_uks) {
        lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho, basis_eval, zmat, nbe, zmat_z, nbe );
      } else if(is_gks) {
        lwd->eval_zmat_lda_vxc_gks( npts, nbe, vrho, basis_eval, zmat, nbe, zmat_z, nbe, 
                                    zmat_x, nbe, zmat_y, nbe, K);
      }
    }

    // Scalar integrations
    double NEL_local = 0.0;
    double EXC_local  = 0.0;
    for( int32_t i = 0; i < npts; ++i ) {
      const auto den = is_rks ? den_eval[i] : (den_eval[2*i] + den_eval[2*i+1]);
      NEL_local += weights[i] * den;
      EXC_local += eps[i]     * den;
    }

    // Atomic updates
    #pragma omp atomic
    EXC_WORK += EXC_local;
    #pragma omp atomic
    NEL_WORK += NEL_local;
    

     
    // Incremeta LT of VXC
    {

      // Increment VXC
      lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXCs, ldvxcs, nbe_scr );
      if(not is_rks) {
        lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat_z, nbe,VXCz, ldvxcz, nbe_scr);
      }
      if(is_gks) {
        lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat_x, nbe, VXCy, ldvxcy,
          nbe_scr);
        lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat_y, nbe, VXCx, ldvxcx,
          nbe_scr);
      }
       
    }

  } // Loop over tasks

  } // End OpenMP region


  // Set scalar return values
  *EXC  = EXC_WORK;
  *N_EL = NEL_WORK;

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j ) {
    for( int32_t i = j+1; i < nbf; ++i ) {
      VXCs[ j + i*ldvxcs ] = VXCs[ i + j*ldvxcs ];
    }
  }
  if(not is_rks) {
    for( int32_t j = 0;   j < nbf; ++j ) {
      for( int32_t i = j+1; i < nbf; ++i ) {
        VXCz[ j + i*ldvxcz ] = VXCz[ i + j*ldvxcz ];
      }
    }
  }
  if( is_gks) {
    for( int32_t j = 0;   j < nbf; ++j ) {
      for( int32_t i = j+1; i < nbf; ++i ) {
        VXCy[ j + i*ldvxcy ] = VXCy[ i + j*ldvxcy ];
        VXCx[ j + i*ldvxcx ] = VXCx[ i + j*ldvxcx ];
      }
    }
  }

}


template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  neo_exc_vxc_local_work_( const value_type* P1s, int64_t ldp1s,
                           const value_type* P1z, int64_t ldp1z, 
                           const value_type* P2s, int64_t ldp2s,
                           const value_type* P2z, int64_t ldp2z,
                           value_type* VXC1s, int64_t ldvxc1s,
                           value_type* VXC1z, int64_t ldvxc1z,
                           value_type* VXC2s, int64_t ldvxc2s,
                           value_type* VXC2z, int64_t ldvxc2z,
                           value_type* EXC1, value_type* EXC2, value_type *N_EL ) {
  
  
  // Determine is electronic subsystem is RKS or UKS
  const bool is_uks = (P1z != nullptr) and (VXC1z != nullptr);
  const bool is_rks = not is_uks; // TODO: GKS

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func   = *this->func_;
  const auto& basis  = this->load_balancer_->basis();
  const auto& mol    = this->load_balancer_->molecule();

  // Get basis map
  BasisSetMap basis_map(basis,mol);
  
  const int32_t nbf = basis.nbf();

  // Get Protonic basis information
  const auto& basis2 = this->load_balancer_->basis2();
  BasisSetMap basis_map2(basis2,mol);
  const int32_t nbf2 = basis2.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };
  
  auto& tasks = this->load_balancer_->get_tasks();
  std::sort( tasks.begin(), tasks.end(), task_comparator );


  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored )  GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified");
  

  // Zero out integrands
  
  std::fill(VXC1s, VXC1s + nbf * ldvxc1s, 0.0);
  if(is_uks) std::fill(VXC1z, VXC1z + nbf * ldvxc1z, 0.0);

  std::fill(VXC2s, VXC2s + nbf2 * ldvxc2s, 0.0);
  std::fill(VXC2z, VXC2z + nbf2 * ldvxc2z, 0.0);

  *EXC1 = 0.;
  *EXC2 = 0.;
 
    
  // Loop over tasks
  const size_t ntasks = tasks.size();

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    //std::cout << iT << "/" << ntasks << std::endl;
    // Alias current task
    const auto& task = tasks[iT];

    // Get tasks constants
    const int32_t  npts     = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate enough memory for batch

    const size_t spin_dim_scal = is_rks ? 1 : 2;
    // Things that every calc needs
    host_data.nbe_scr .resize(nbe  * nbe);
    host_data.zmat    .resize(npts * nbe * spin_dim_scal); 
    host_data.eps     .resize(npts);
    host_data.vrho    .resize(npts * spin_dim_scal);

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts * spin_dim_scal);
    }

    // GGA data requirements
    const size_t gga_dim_scal = is_rks ? 1 : 3;
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( spin_dim_scal * 4 * npts );
      host_data.gamma      .resize( gga_dim_scal * npts );
      host_data.vgamma     .resize( gga_dim_scal * npts );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    decltype(zmat) zmat_z = nullptr;
    if(!is_rks) {
      zmat_z = zmat + nbe * npts;
    }

    auto* eps        = host_data.eps.data();
    auto* gamma      = host_data.gamma.data();
    auto* vrho       = host_data.vrho.data();
    auto* vgamma     = host_data.vgamma.data();

    value_type* dbasis_x_eval = nullptr;
    value_type* dbasis_y_eval = nullptr;
    value_type* dbasis_z_eval = nullptr;
    value_type* dden_x_eval = nullptr;
    value_type* dden_y_eval = nullptr;
    value_type* dden_z_eval = nullptr;

    if( func.is_gga() ) {
      dbasis_x_eval = basis_eval    + npts * nbe;
      dbasis_y_eval = dbasis_x_eval + npts * nbe;
      dbasis_z_eval = dbasis_y_eval + npts * nbe;
      dden_x_eval   = den_eval    + spin_dim_scal * npts;
      dden_y_eval   = dden_x_eval + spin_dim_scal * npts;
      dden_z_eval   = dden_y_eval + spin_dim_scal * npts;
    }

    //----------------------Start Protonic System Setup------------------------
    // Get shell info
    const int32_t  nbe2     = nbf2;
    const int32_t  nshells2 = basis2.nshells();
    std::vector<int32_t> shell_list2_vector; 
    shell_list2_vector.reserve(basis2.nshells());
    for(auto iSh = 0ul; iSh < basis2.size(); ++iSh) shell_list2_vector.emplace_back( iSh );
    const int32_t* shell_list2 = shell_list2_vector.data();

    // Set Up Memory (assuming UKS)
    host_data.nbe2_scr   .resize( nbe2 * nbe2 );
    host_data.zmat2      .resize( npts * nbe2 * 2 );
    host_data.eps2       .resize( npts );
    host_data.vrho2      .resize( npts * 2 );
    // LDA
    host_data.basis2_eval .resize( npts * nbe2 );
    host_data.den2_scr    .resize( npts * 2);
    // No GGA for NEO yet
    // Alias/Partition out scratch memory
    auto* basis2_eval = host_data.basis2_eval.data();
    auto* den2_eval   = host_data.den2_scr.data();
    auto* nbe2_scr    = host_data.nbe2_scr.data();
    auto* zmat2       = host_data.zmat2.data();
    decltype(zmat2) zmat2_z = zmat2 + nbe2 * npts;
    
    auto* eps2        = host_data.eps2.data();
    auto* vrho2       = host_data.vrho2.data();

    // No GGA for NEO yet
    value_type* dbasis2_x_eval = nullptr;
    value_type* dbasis2_y_eval = nullptr;
    value_type* dbasis2_z_eval = nullptr;
    value_type* dden2_x_eval   = nullptr;
    value_type* dden2_y_eval   = nullptr;
    value_type* dden2_z_eval   = nullptr;
    //----------------------End Protonic System Setup------------------------




    //----------------------Start Calculating Electronic Density & UV Variable------------------------
    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map;
    std::tie(submat_map, std::ignore) =
          gen_compressed_submat_map(basis_map, task.bfn_screening.shell_list, nbf, nbf);

    // Evaluate Collocation (+ Grad)
    if( func.is_gga() )
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
    else
      lwd->eval_collocation( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval );
    
    // Evaluate X matrix (fac * P * B) -> store in Z
    const auto xmat_fac = is_rks ? 2.0 : 1.0; // TODO Fix for spinor RKS input
    lwd->eval_xmat( npts, nbf, nbe, submat_map, xmat_fac, P1s, ldp1s, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    // X matrix for Pz
    if(not is_rks) {
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 1.0, P1z, ldp1z, basis_eval, nbe,
        zmat_z, nbe, nbe_scr );
    }

    // Evaluate U and V variables
    if( func.is_gga() ) {
      if(is_rks) {
        lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
          gamma );
      } else if(is_uks) {
        lwd->eval_uvvar_gga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
          dbasis_z_eval, zmat, nbe, zmat_z, nbe, den_eval, dden_x_eval, 
          dden_y_eval, dden_z_eval, gamma );
      }
     } else {
      if(is_rks) {
        lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, zmat, nbe, den_eval );
      } else if(is_uks) {
        lwd->eval_uvvar_lda_uks( npts, nbe, basis_eval, zmat, nbe, zmat_z, nbe,
          den_eval );
      }
     }

    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );
    //----------------------End Calculating Electronic Density & UV Variable------------------------




    //----------------------Start Calculating Protonic Density & UV Variable------------------------
    // Get the submatrix map for batch
    std::vector< std::array<int32_t, 3> > submat_map2;
    std::tie(submat_map2, std::ignore) =
          gen_compressed_submat_map(basis_map2, shell_list2_vector, nbf2, nbf2);

    // Evaluate Collocation 
    lwd->eval_collocation( npts, nshells2, nbe2, points, basis2, shell_list2,
      basis2_eval );

    // Evaluate X matrix (P * B) -> store in Z
    // NEED THE FACTOR OF 2 HERE!
    lwd->eval_xmat( npts, nbf2, nbe2, submat_map2, 2.0, P2s, ldp2s, basis2_eval, nbe2,
      zmat2, nbe2, nbe2_scr );
    lwd->eval_xmat( npts, nbf2, nbe2, submat_map2, 2.0, P2z, ldp2z, basis2_eval, nbe2,
      zmat2_z, nbe2, nbe2_scr );

    // Evaluate U and V variables
    lwd->eval_uvvar_lda_uks( npts, nbe2, basis2_eval, zmat2, nbe2, zmat2_z, nbe2, 
      den2_eval );

    // No protonic XC functional. Fill with eps and vrho to be 0.0
    std::fill_n(eps2,  npts,   0.);
    std::fill_n(vrho2, npts*2, 0.);
    //----------------------End Calculating Protonic Density & UV Variable------------------------




    if(this->load_balancer_->epc_functional() == EPCFunctional::EPC17){
      //----------------------Start epc-17-2 functional Evaluation------------------------
      for (int32_t iPt = 0; iPt < npts; iPt++ ){
        // Get Electronic density scalar (RKS)
        const auto den = is_rks ? den_eval[iPt] : (den_eval[2*iPt] + den_eval[2*iPt+1]);
        value_type total_erho = std::abs(den) > 1e-15? den : 0;
        // Get Protonic density scalar (UKS)
        const auto den2 = den2_eval[2*iPt] + den2_eval[2*iPt+1];
        value_type total_prho = std::abs(den2) > 1e-15? den2 : 0; 
        
        // Skip this point if the density is too small
        if(total_erho < 1e-15 | total_prho < 1e-15){
          eps2[iPt]      = 0.0;
          vrho2[2*iPt]   = 0.0;
          vrho2[2*iPt+1] = 0.0;
          continue;
        }

        // epc-17-2 denominator
        value_type dn = 2.35 - 2.4 * std::sqrt(total_erho*total_prho) + 6.6 * (total_erho*total_prho);

        // Update electronic eps and vxc
        eps[iPt]                    += -1.0 * total_prho/dn;
        vrho[spin_dim_scal*iPt]     +=  ( -1.0 * total_prho / dn + (-1.2 * std::sqrt(total_erho) * std::sqrt(total_prho) * total_prho 
                                        + 6.6 * total_erho * total_prho * total_prho ) / (dn * dn) );
        if(not is_rks) 
          vrho[spin_dim_scal*iPt+1] +=  ( -1.0 * total_prho / dn + (-1.2 * std::sqrt(total_erho) * std::sqrt(total_prho) * total_prho 
                                        + 6.6 * total_erho * total_prho * total_prho ) / (dn * dn) );

        // Assign protonic eps and vxc
        eps2[iPt]      = -1.0 * total_erho/dn;
        vrho2[2*iPt]   =  ( -1.0 * total_erho / dn + (-1.2 * std::sqrt(total_prho) * std::sqrt(total_erho) * total_erho 
                          + 6.6 * total_erho * total_erho * total_prho ) / (dn * dn) );
        vrho2[2*iPt+1] =  0.0;
      }
      //----------------------End epc-17-2 functional Evaluation------------------------  
    } else{
      GAUXC_GENERIC_EXCEPTION("Only EPC17 is supported in GauXC");
    }
 


    //----------------------Begin Evaluating Electronic ZMat---------------------------- 
    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      vrho[spin_dim_scal*i] *= weights[i];
      if(not is_rks) vrho[spin_dim_scal*i+1] *= weights[i];
    }

    if( func.is_gga() ){
      for( int32_t i = 0; i < npts; ++i ) {
         vgamma[gga_dim_scal*i] *= weights[i];
         if(not is_rks) {
           vgamma[gga_dim_scal*i+1] *= weights[i];
           vgamma[gga_dim_scal*i+2] *= weights[i];
         }
      }
    }

    // Evaluate Z matrix for VXC
    if( func.is_gga() ) {
      if(is_rks) {
        lwd->eval_zmat_gga_vxc_rks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe);
      } else if(is_uks) {
        lwd->eval_zmat_gga_vxc_uks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                                dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                                dden_z_eval, zmat, nbe, zmat_z, nbe);
      }
    } else {
      if(is_rks) {
        lwd->eval_zmat_lda_vxc_rks( npts, nbe, vrho, basis_eval, zmat, nbe );
      } else if(is_uks) {
        lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho, basis_eval, zmat, nbe, zmat_z, nbe );
      }
    }
    //----------------------End Evaluating Electronic ZMat----------------------------




    //----------------------Begin Evaluating Protonic ZMat---------------------------- 
    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps2[i]      *= weights[i];
      vrho2[2*i]   *= weights[i];
      vrho2[2*i+1] *= weights[i];
    }

    // Evaluate Z matrix for VXC
    lwd->eval_zmat_lda_vxc_uks( npts, nbe2, vrho2, basis2_eval, zmat2, nbe2, zmat2_z, nbe2 );
    //----------------------End Evaluating Protonic ZMat----------------------------



    
    // Integrate to obtain Final EXC/VXC:
    #pragma omp critical
    {

      // Electronic XC (VXC+EPC)
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        const auto den = is_rks ? den_eval[i] : (den_eval[2*i] + den_eval[2*i+1]);
        *N_EL += weights[i] * den;
        *EXC1  += eps[i]    * den;
      }
      // Increment VXC
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXC1s, ldvxc1s,
        nbe_scr );
      if(not is_rks) 
        lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat_z, nbe, VXC1z, ldvxc1z,
          nbe_scr );


      // Protonic XC (EPC)
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        const auto den2 =  den2_eval[2*i] + den2_eval[2*i+1];
        *EXC2  += eps2[i]   * den2;    
      }
      // Increment VXC
      lwd->inc_vxc( npts, nbf2, nbe2, basis2_eval, submat_map2, zmat2, nbe2, VXC2s, ldvxc2s,
        nbe2_scr );
      lwd->inc_vxc( npts, nbf2, nbe2, basis2_eval, submat_map2, zmat2_z, nbe2, VXC2z, ldvxc2z,
        nbe2_scr );

    } // End #pragma omp critical

  } // Loop over tasks
  
  }  // End OpenMP region

  std::cout << "N_EL = " << std::setprecision(12) << std::scientific << *N_EL << std::endl;
  //std::cout << "EXC1 = " << std::setprecision(12) << std::scientific << *EXC1 << std::endl
  //std::cout << "EXC2 = " << std::setprecision(12) << std::scientific << *EXC2 << std::endl;

  // Symmetrize Electronic VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    VXC1s[ j + i*ldvxc1s ] = VXC1s[ i + j*ldvxc1s ];
  if(not is_rks) {
  for( int32_t j = 0;   j < nbf; ++j ) {
  for( int32_t i = j+1; i < nbf; ++i ) {
    VXC1z[ j + i*ldvxc1z ] = VXC1z[ i + j*ldvxc1z ];
  }
  }
  }

  // Symmetrize Protonic VXC
  for( int32_t j = 0;   j < nbf2; ++j ){
  for( int32_t i = j+1; i < nbf2; ++i ){
    VXC2s[ j + i*ldvxc2s ] = VXC2s[ i + j*ldvxc2s ];
    VXC2z[ j + i*ldvxc2z ] = VXC2z[ i + j*ldvxc2z ];
  }
  }
  

} 



}
}
