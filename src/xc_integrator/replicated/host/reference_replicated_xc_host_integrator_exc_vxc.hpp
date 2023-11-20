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
#include "host/blas.hpp"
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
  this->load_balancer_->get_tasks();

  // Temporary electron count to judge integrator accuracy
  value_type N_EL;

  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    //exc_vxc_local_work_( P, ldp, VXC, ldvxc, EXC, &N_EL );
    exc_vxc_local_work_( P, ldp, nullptr, 0, VXC, ldvxc, nullptr, 0, EXC, &N_EL );
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
                      value_type* EXC ) {

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
    exc_vxc_local_work_( Ps, ldps, Pz, ldpz, VXCs, ldvxcs, VXCz, ldvxcz, EXC, &N_EL );
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
  exc_vxc_local_work_( const value_type* Ps, int64_t ldps,
                       const value_type* Pz, int64_t ldpz,
                       value_type* VXCs, int64_t ldvxcs,
                       value_type* VXCz, int64_t ldvxcz, 
                       value_type* EXC, value_type *N_EL ) {


  const bool is_uks = (Pz != nullptr) and (VXCz != nullptr);
  const bool is_rks = not is_uks; // TODO: GKS

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();

  const bool needs_laplacian = func.is_mgga() ? true : false; // TODO: Needs Laplacian Check
  
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
  if(is_uks) {
    for( auto j = 0; j < nbf; ++j ) {
      for( auto i = 0; i < nbf; ++i ) {
        VXCz[i + j*ldvxcz] = 0.;
      }
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

    const size_t spin_dim_scal = is_rks ? 1 : 2; 
    const size_t mgga_dim_scal = func.is_mgga() ? 4 : 1; // basis + d1basis

    // Things that every calc needs
    host_data.nbe_scr .resize(nbe  * nbe);
    host_data.zmat    .resize(npts * nbe * spin_dim_scal * mgga_dim_scal); 
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
      // TODO: Add check for Laplacian-dependent functionals
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
    if(!is_rks) {
      zmat_z = zmat + mgga_dim_scal * nbe * npts;
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
    if( func.is_mgga() )
      func.eval_exc_vxc( npts, den_eval, gamma, lapl, tau, eps, vrho, vgamma, vlapl, vtau);
    else if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );

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
      }
    } else {
      if(is_rks) {
        lwd->eval_zmat_lda_vxc_rks( npts, nbe, vrho, basis_eval, zmat, nbe );
      } else if(is_uks) {
        lwd->eval_zmat_lda_vxc_uks( npts, nbe, vrho, basis_eval, zmat, nbe, zmat_z, nbe );
      }
    }


    // Incremeta LT of VXC
    #pragma omp critical
    {
      // Scalar integrations
      for( int32_t i = 0; i < npts; ++i ) {
        const auto den = is_rks ? den_eval[i] : (den_eval[2*i] + den_eval[2*i+1]);
        *N_EL += weights[i] * den;
        *EXC  += eps[i]     * den;
      }

      // Increment VXC
      lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat, nbe, VXCs, ldvxcs, nbe_scr );
      if(not is_rks) {
        lwd->inc_vxc( mgga_dim_scal * npts, nbf, nbe, basis_eval, submat_map, zmat_z, nbe,VXCz, ldvxcz, nbe_scr);
      }

    }

  } // Loop over tasks

  } // End OpenMP region

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

} 

}
}
