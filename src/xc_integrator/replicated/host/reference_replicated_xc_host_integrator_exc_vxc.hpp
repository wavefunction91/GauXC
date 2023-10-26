/**
 * GauXC Copyright (c) 2020-2023, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 *
 * See LICENSE.txt for details
 */
#pragma once

#include "reference_replicated_xc_host_integrator.hpp"
#include "integrator_util/integrator_common.hpp"
#include "host/local_host_work_driver.hpp"
#include <stdexcept>

namespace GauXC  {
namespace detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* VXC, int64_t ldvxc,
                 value_type* EXC ) {

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
  auto& tasks = this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( basis, P, ldp, VXC, ldvxc, EXC, &N_EL,
      tasks.begin(), tasks.end() );
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
  eval_exc_vxc_( int64_t m, int64_t n, const value_type* Pscalar,
                      int64_t ldpscalar,
                      const value_type* Pz,
                      int64_t ldpz,
                      value_type* VXCscalar, int64_t ldvxcscalar,
                      value_type* VXCz, int64_t ldvxcz,
                      value_type* EXC ) {

  const auto& basis = this->load_balancer_->basis();

  // Check that P / VXC are sane
  const int64_t nbf = basis.nbf();
  if( m != n )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Be Square");
  if( m != nbf )
    GAUXC_GENERIC_EXCEPTION("P/VXC Must Have Same Dimension as Basis");
  if( ldpscalar < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPSCALAR");
  if( ldpz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDPZ");
  if( ldvxcscalar < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCSCALAR");
  if( ldvxcz < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDVXCZ");

  // Get Tasks
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_vxc_local_work_( Pscalar, ldpscalar, Pz, ldpz,  VXCscalar, ldvxcscalar, VXCz, ldvxcz, EXC, &N_EL );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    this->reduction_driver_->allreduce_inplace( VXCscalar, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( VXCz, nbf*nbf, ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( EXC,   1    , ReductionOp::Sum );
    this->reduction_driver_->allreduce_inplace( &N_EL, 1    , ReductionOp::Sum );

  });


}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_vxc_local_work_(const basis_type& basis,  const value_type* P, int64_t ldp, 
    value_type* VXC, int64_t ldvxc, value_type* EXC, 
    value_type* N_EL, task_iterator task_begin, task_iterator task_end ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& mol   = this->load_balancer_->molecule();

  // Get basis map
  BasisSetMap basis_map(basis,mol);

  const int32_t nbf = basis.nbf();

  // Sort tasks on size (XXX: maybe doesnt matter?)
  auto task_comparator = []( const XCTask& a, const XCTask& b ) {
    return (a.points.size() * a.bfn_screening.nbe) > (b.points.size() * b.bfn_screening.nbe);
  };

  std::sort( task_begin, task_end, task_comparator );


  // Check that Partition Weights have been calculated
  auto& lb_state = this->load_balancer_->state();
  if( not lb_state.modified_weights_are_stored ) {
    GAUXC_GENERIC_EXCEPTION("Weights Have Not Beed Modified"); 
  }

  // Zero out integrands
  for( auto j = 0; j < nbf; ++j )
  for( auto i = 0; i < nbf; ++i ) 
    VXC[i + j*ldvxc] = 0.;
  *EXC = 0.;


  // Loop over tasks
  const size_t ntasks = std::distance(task_begin, task_end);

  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

    //std::cout << iT << "/" << ntasks << std::endl;
    // Alias current task
    const auto& task = *(task_begin + iT);

    // Get tasks constants
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate enough memory for batch

    // Things that every calc needs
    host_data.nbe_scr .resize( nbe * nbe  );
    host_data.zmat    .resize( npts * nbe );
    host_data.eps     .resize( npts );
    host_data.vrho    .resize( npts );

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts );
    }

    // GGA data requirements
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( 4 * npts );
      host_data.gamma      .resize( npts );
      host_data.vgamma     .resize( npts );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

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
      dden_x_eval   = den_eval    + npts;
      dden_y_eval   = dden_x_eval + npts;
      dden_z_eval   = dden_y_eval + npts;
    }


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


    // Evaluate X matrix (P * B) -> store in Z
    lwd->eval_xmat( npts, nbf, nbe, submat_map, P, ldp, basis_eval, nbe,
      zmat, nbe, nbe_scr );


    // Evaluate U and V variables
    if( func.is_gga() )
      lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
        dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
        gamma );
     else
      lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, zmat, nbe, den_eval );

    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      vrho[i] *= weights[i];
    }

    if( func.is_gga() )
      for( int32_t i = 0; i < npts; ++i ) vgamma[i] *= weights[i];




    // Evaluate Z matrix for VXC
    if( func.is_gga() )
      lwd->eval_zmat_gga_vxc_rks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                              dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                              dden_z_eval, zmat, nbe); 
    else
      lwd->eval_zmat_lda_vxc_rks( npts, nbe, vrho, basis_eval, zmat, nbe ); 


    // Incremeta LT of VXC
    #pragma omp critical
    {
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        *N_EL += weights[i] * den_eval[i];
        *EXC  += eps[i]     * den_eval[i];
      }

      // Increment VXC
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXC, ldvxc,
        nbe_scr );
    }

  } // Loop over tasks 

  } // End OpenMP region

  //std::cout << "N_EL = " << std::setprecision(12) << std::scientific << *N_EL << std::endl;

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j )
  for( int32_t i = j+1; i < nbf; ++i )
    VXC[ j + i*nbf ] = VXC[ i + j*nbf ];

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_vxc_local_work_( const value_type* Pscalar, int64_t ldpscalar,
                            const value_type* Pz, int64_t ldpz,
                            value_type* VXCscalar, int64_t ldvxcscalar,
                            value_type* VXCz, int64_t ldvxcz, value_type* EXC, value_type *N_EL ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();

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
    VXCscalar[i + j*ldvxcscalar] = 0.;
    VXCz[i + j*ldvxcz] = 0.;
  }
  }
  *EXC = 0.;
 
    
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
    const int32_t  npts    = task.points.size();
    const int32_t  nbe     = task.bfn_screening.nbe;
    const int32_t  nshells = task.bfn_screening.shell_list.size();

    const auto* points      = task.points.data()->data();
    const auto* weights     = task.weights.data();
    const int32_t* shell_list = task.bfn_screening.shell_list.data();

    // Allocate enough memory for batch

    // Things that every calc needs
    host_data.nbe_scr .resize( nbe * nbe * 2 );
    host_data.zmat    .resize( npts * nbe * 2);
    host_data.eps     .resize( npts );
    host_data.vrho    .resize( npts * 2);

    // LDA data requirements
    if( func.is_lda() ){
      host_data.basis_eval .resize( npts * nbe );
      host_data.den_scr    .resize( npts * 2);
    }

    // GGA data requirements
    if( func.is_gga() ){
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.den_scr    .resize( 2 * 4 * npts );
      host_data.gamma      .resize( 3 * npts );
      host_data.vgamma     .resize( 3 * npts );
    }

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

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
      dden_x_eval   = den_eval    + 2*npts;
      dden_y_eval   = dden_x_eval + 2*npts;
      dden_z_eval   = dden_y_eval + 2*npts;
    }


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


    // Evaluate X matrix (P * B) -> store in Z
    lwd->eval_xmat( npts, nbf, nbe, submat_map, Pscalar, ldpscalar, basis_eval, nbe,
      zmat, nbe, nbe_scr );

    lwd->eval_xmat( npts, nbf, nbe, submat_map, Pz, ldpz, basis_eval, nbe,
      zmat + npts*nbe, nbe, nbe_scr + nbe * nbe);


    // Evaluate U and V variables
    if( func.is_gga() )
      lwd->eval_uvvar_gga_uks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
        dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
        gamma );
     else
      lwd->eval_uvvar_lda_uks( npts, nbe, basis_eval, zmat, nbe, den_eval );
    
    // Evaluate XC functional
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

    // Factor weights into XC results
    for( int32_t i = 0; i < npts; ++i ) {
      eps[i]  *= weights[i];
      vrho[2*i] *= weights[i];
      vrho[2*i+1] *= weights[i];
    }

    if( func.is_gga() ){
      for( int32_t i = 0; i < npts; ++i ) {
         vgamma[3*i] *= weights[i];
         vgamma[3*i+1] *= weights[i];
         vgamma[3*i+2] *= weights[i];
      }
    }



    // Evaluate Z matrix for VXC
    if( func.is_gga() )
      lwd->eval_zmat_gga_vxc_uks( npts, nbe, vrho, vgamma, basis_eval, dbasis_x_eval,
                              dbasis_y_eval, dbasis_z_eval, dden_x_eval, dden_y_eval,
                              dden_z_eval, zmat, nbe);
    else
      lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho, basis_eval, zmat, nbe );


    // Incremeta LT of VXC
    #pragma omp critical
    {
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        *N_EL += weights[i] * (den_eval[2*i] +  den_eval[2*i+1]);
        *EXC  += eps[i]     * (den_eval[2*i] +  den_eval[2*i+1]);
      }

      // Increment VXC
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXCscalar, ldvxcscalar,
        nbe_scr );
      lwd->inc_vxc( npts, nbf, nbe, basis_eval, submat_map, zmat+ npts*nbe, nbe, VXCz, ldvxcz,
        nbe_scr + nbe * nbe);

    }

  } // Loop over tasks

  } // End OpenMP region

  //std::cout << "N_EL = " << std::setprecision(12) << std::scientific << *N_EL << std::endl;

  // Symmetrize VXC
  for( int32_t j = 0;   j < nbf; ++j ) {
  for( int32_t i = j+1; i < nbf; ++i ) {
    VXCscalar[ j + i*nbf ] = VXCscalar[ i + j*nbf ];
    VXCz[ j + i*nbf ] = VXCz[ i + j*nbf ];
  }
  }

} 

}
}
