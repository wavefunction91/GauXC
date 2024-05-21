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

namespace GauXC::detail {

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  eval_exc_grad_( int64_t m, int64_t n, const value_type* P,
                 int64_t ldp, value_type* EXC_GRAD ) { 
                 
                 
  const auto& basis = this->load_balancer_->basis();

  // Check that P is sane
  const int64_t nbf = basis.nbf();
  if( m != n ) 
    GAUXC_GENERIC_EXCEPTION("P Must Be Square");
  if( m != nbf ) 
    GAUXC_GENERIC_EXCEPTION("P Must Have Same Dimension as Basis");
  if( ldp < nbf )
    GAUXC_GENERIC_EXCEPTION("Invalid LDP");
                 
                 
  // Get Tasks
  this->load_balancer_->get_tasks();
                 
  // Compute Local contributions to EXC / VXC
  this->timer_.time_op("XCIntegrator.LocalWork", [&](){
    exc_grad_local_work_( P, ldp, EXC_GRAD );
  });


  // Reduce Results
  this->timer_.time_op("XCIntegrator.Allreduce", [&](){

    if( not this->reduction_driver_->takes_host_memory() )
      GAUXC_GENERIC_EXCEPTION("This Module Only Works With Host Reductions");

    const int natoms = this->load_balancer_->molecule().natoms();
    this->reduction_driver_->allreduce_inplace( EXC_GRAD, 3*natoms, ReductionOp::Sum );
  });

}

template <typename ValueType>
void ReferenceReplicatedXCHostIntegrator<ValueType>::
  exc_grad_local_work_( const value_type* P, int64_t ldp, value_type* EXC_GRAD ) {

  // Cast LWD to LocalHostWorkDriver
  auto* lwd = dynamic_cast<LocalHostWorkDriver*>(this->local_work_driver_.get());

  // Setup Aliases
  const auto& func  = *this->func_;
  const auto& basis = this->load_balancer_->basis();
  const auto& mol   = this->load_balancer_->molecule();

  // MGGA constants
  const size_t mmga_dim_scal = func.is_mgga() ? 4 : 1;
  const bool needs_laplacian = func.is_mgga() ? true : false; // TODO: Check for Laplacian dependence
							      //
  // Get basis map
  BasisSetMap basis_map(basis,mol);

  const int32_t nbf = basis.nbf();
  const int32_t natoms = mol.natoms();

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
  for( auto i = 0; i < 3*natoms; ++i ) {
    EXC_GRAD[i] = 0.;
  }

  // Loop over tasks
  const size_t ntasks = tasks.size();
  #pragma omp parallel
  {

  XCHostData<value_type> host_data; // Thread local host data

  #pragma omp for schedule(dynamic)
  for( size_t iT = 0; iT < ntasks; ++iT ) {

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
    host_data.nbe_scr .resize( nbe * nbe  );
    host_data.eps     .resize( npts );
    host_data.vrho    .resize( npts );
    host_data.den_scr .resize( 4 * npts );

    if( func.is_lda() ) {
      host_data.basis_eval .resize( 4 * npts * nbe );
      host_data.zmat       .resize( npts * nbe );
    }

    if( func.is_gga() ){
      host_data.basis_eval .resize( 10 * npts * nbe );
      host_data.zmat       .resize( 4  * npts * nbe );
      host_data.gamma      .resize( npts );
      host_data.vgamma     .resize( npts );
    }

#if 0
    if( func.is_mgga() ) {
      host_data.basis_eval .resize( 11 * npts * nbe ); // basis + grad(3) + hess(6) + lapl
      host_data.zmat       .resize(  7 * npts * nbe ); // basis + grad(3) + grad(3)
      host_data.mmat       .resize( npts * nbe );
      host_data.gamma      .resize( npts );
      host_data.vgamma     .resize( npts );
      host_data.tau        .resize( npts );
      host_data.vtau       .resize( npts );
      if ( needs_laplacian ) {
	host_data.basis_eval.resize( 24 * npts * nbe );
	host_data.lapl      .resize( npts );
	host_data.vlapl     .resize( npts );
      }
    }
#endif

    // Alias/Partition out scratch memory
    auto* basis_eval = host_data.basis_eval.data();
    auto* den_eval   = host_data.den_scr.data();
    auto* nbe_scr    = host_data.nbe_scr.data();
    auto* zmat       = host_data.zmat.data();

    auto* zmat_x = zmat   + npts*nbe;
    auto* zmat_y = zmat_x + npts*nbe;
    auto* zmat_z = zmat_y + npts*nbe;

    auto* eps        = host_data.eps.data();
    auto* gamma      = host_data.gamma.data();
    auto* vrho       = host_data.vrho.data();
    auto* vgamma     = host_data.vgamma.data();

#if 0
    auto* tau        = host_data.tau.data();
    auto* lapl       = host_data.lapl.data();
    auto* vtau       = host_data.vtau.data();
    auto* vlapl      = host_data.vlapl.data();
    auto* mmat_x      = mmat;
    auto* mmat_y      = mmat_x + npts * nbe;
    auto* mmat_z      = mmat_y + npts * nbe;
#endif

    auto* dbasis_x_eval = basis_eval    + npts * nbe;
    auto* dbasis_y_eval = dbasis_x_eval + npts * nbe;
    auto* dbasis_z_eval = dbasis_y_eval + npts * nbe;
    auto* dden_x_eval   = den_eval    + npts;
    auto* dden_y_eval   = dden_x_eval + npts;
    auto* dden_z_eval   = dden_y_eval + npts;

    value_type* d2basis_xx_eval = nullptr;
    value_type* d2basis_xy_eval = nullptr;
    value_type* d2basis_xz_eval = nullptr;
    value_type* d2basis_yy_eval = nullptr;
    value_type* d2basis_yz_eval = nullptr;
    value_type* d2basis_zz_eval = nullptr;
#if 0
    value_type* lbasis_eval = nullptr;
    value_type* d3basis_xxx_eval = nullptr;
    value_type* d3basis_xxy_eval = nullptr;
    value_type* d3basis_xxz_eval = nullptr;
    value_type* d3basis_xyy_eval = nullptr;
    value_type* d3basis_xyz_eval = nullptr;
    value_type* d3basis_xzz_eval = nullptr;
    value_type* d3basis_yyy_eval = nullptr;
    value_type* d3basis_yyz_eval = nullptr;
    value_type* d3basis_yzz_eval = nullptr;
    value_type* d3basis_zzz_eval = nullptr;
    value_type* dlbasis_x_eval = nullptr;
    value_type* dlbasis_y_eval = nullptr;
    value_type* dlbasis_z_eval = nullptr;
#endif

    if( func.is_gga() ) {
      d2basis_xx_eval = dbasis_z_eval   + npts * nbe;
      d2basis_xy_eval = d2basis_xx_eval + npts * nbe;
      d2basis_xz_eval = d2basis_xy_eval + npts * nbe;
      d2basis_yy_eval = d2basis_xz_eval + npts * nbe;
      d2basis_yz_eval = d2basis_yy_eval + npts * nbe;
      d2basis_zz_eval = d2basis_yz_eval + npts * nbe;
    }

#if 0
    if( func.is_mgga() ) {
      d2basis_xx_eval = dbasis_z_eval   + npts * nbe;
      d2basis_xy_eval = d2basis_xx_eval + npts * nbe;
      d2basis_xz_eval = d2basis_xy_eval + npts * nbe;
      d2basis_yy_eval = d2basis_xz_eval + npts * nbe;
      d2basis_yz_eval = d2basis_yy_eval + npts * nbe;
      d2basis_zz_eval = d2basis_yz_eval + npts * nbe;
      if ( true ) {
	lbasis_eval   = d2basis_zz_eval + npts * nbe;
	d3basis_xxx_eval = lbasis_eval + npts * nbe;
	d3basis_xxy_eval = d3basis_xxx_eval + npts * nbe;
	d3basis_xxz_eval = d3basis_xxy_eval + npts * nbe;
	d3basis_xyy_eval = d3basis_xxz_eval + npts * nbe;
	d3basis_xyz_eval = d3basis_xyy_eval + npts * nbe;
	d3basis_xzz_eval = d3basis_xyz_eval + npts * nbe;
	d3basis_yyy_eval = d3basis_xzz_eval + npts * nbe;
	d3basis_yyz_eval = d3basis_yyy_eval + npts * nbe;
	d3basis_yzz_eval = d3basis_yyz_eval + npts * nbe;
	d3basis_zzz_eval = d3basis_yzz_eval + npts * nbe;
        dlbasis_x_eval = d3basis_zzz_eval + npts * nbe;
	dlbasis_y_eval = dlbasis_x_eval + npts * nbe;
	dlbasis_z_eval = dlbasis_y_eval + npts * nbe;
      }
    }
#endif


    // Get the submatrix map for batch
    auto [submat_map, foo] = 
      gen_compressed_submat_map( basis_map, task.bfn_screening.shell_list, nbf, nbf );

    // Evaluate Collocation Gradient (+ Hessian)
#if 0
    if( func.is_mgga() ) {
      lwd->eval_collocation_der3( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval,
        d2basis_xy_eval, d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval,
        d2basis_zz_eval, d3basis_xxx_eval, d3basis_xxy_eval, d3basis_xxz_eval,
	d3basis_xyy_eval, d3basis_xyz_eval, d3basis_xzz_eval, d3basis_yyy_eval,
	d3basis_yyz_eval, d3basis_yzz_eval, d3basis_zzz_eval);

    }
    else if( func.is_gga() )
#endif
    if( func.is_gga() )
      lwd->eval_collocation_hessian( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval,
        d2basis_xy_eval, d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval,
        d2basis_zz_eval );
    else
      lwd->eval_collocation_gradient( npts, nshells, nbe, points, basis, shell_list, 
        basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );


    // Evaluate X matrix (2 * P * B/Bx/By/Bz) -> store in Z
    // XXX: This assumes that bfn + gradients are contiguous in memory
    if( func.is_gga() or func.is_mgga() ) {
      lwd->eval_xmat( 4*npts, nbf, nbe, submat_map, 2.0, P, ldp, basis_eval, nbe,
        zmat, nbe, nbe_scr );
    } else {
      lwd->eval_xmat( npts, nbf, nbe, submat_map, 2.0, P, ldp, basis_eval, nbe,
        zmat, nbe, nbe_scr );
    }

    // Evaluate U and V variables
#if 0
    if( func.is_mgga() ) {
      if ( needs_laplacian ) {
        blas::lacpy( 'A', nbe, npts, d2basis_xx_eval, nbe, lbasis_eval, nbe );
        blas::axpy( nbe * npts, 1., d2basis_yy_eval, 1, lbasis_eval, 1);
        blas::axpy( nbe * npts, 1., d2basis_zz_eval, 1, lbasis_eval, 1);

	blas::lacpy( 'A', nbe, npts, d3basis_xxx_eval, nbe, dlbasis_x_eval, nbe );
        blas::axpy( nbe * npts, 1., d3basis_xyy_eval, 1, dlbasis_x_eval, 1);
        blas::axpy( nbe * npts, 1., d3basis_xzz_eval, 1, dlbasis_x_eval, 1);

	blas::lacpy( 'A', nbe, npts, d3basis_xxy_eval, nbe, dlbasis_y_eval, nbe );
        blas::axpy( nbe * npts, 1., d3basis_yyy_eval, 1, dlbasis_y_eval, 1);
        blas::axpy( nbe * npts, 1., d3basis_yzz_eval, 1, dlbasis_y_eval, 1);

	blas::lacpy( 'A', nbe, npts, d3basis_xxz_eval, nbe, dlbasis_z_eval, nbe );
        blas::axpy( nbe * npts, 1., d3basis_yyz_eval, 1, dlbasis_z_eval, 1);
        blas::axpy( nbe * npts, 1., d3basis_zzz_eval, 1, dlbasis_z_eval, 1);
      }
      lwd->eval_uvvar_mgga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
        dbasis_z_eval, lbasis_eval, zmat, nbe, mmat_x, mmat_y, mmat_z, nbe, 
	den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
        gamma, tau, lapl );
    }
    else if( func.is_gga() )
#endif
    if( func.is_gga() )
      lwd->eval_uvvar_gga_rks( npts, nbe, basis_eval, dbasis_x_eval, dbasis_y_eval,
        dbasis_z_eval, zmat, nbe, den_eval, dden_x_eval, dden_y_eval, dden_z_eval,
        gamma );
     else
      lwd->eval_uvvar_lda_rks( npts, nbe, basis_eval, zmat, nbe, den_eval );


    // Evaluate XC functional
#if 0
    if( func.is_mgga() )
      func.eval_exc_vxc( npts, den_eval, gamma, lapl, tau, eps, vrho, vgamma, vlapl, vtau );
    else if(func.is_gga() )
#endif
    if( func.is_gga() )
      func.eval_exc_vxc( npts, den_eval, gamma, eps, vrho, vgamma );
    else
      func.eval_exc_vxc( npts, den_eval, eps, vrho );


    // Increment EXC Gradient
    size_t bf_off = 0;
    for( auto ish = 0; ish < nshells; ++ish ) {
      const int sh_idx = shell_list[ish];
      const int sh_sz  = basis[sh_idx].size();
      const int iAt    = basis_map.shell_to_center( sh_idx );

      double g_acc_x(0), g_acc_y(0), g_acc_z(0);
      for( int ibf = 0, mu = bf_off; ibf < sh_sz; ++ibf, ++mu )
      for( int ipt = 0; ipt < npts; ++ipt ) {

	      const int32_t mu_i = mu + ipt*nbe;

	      // LDA Contributions
        const double vrho_ipt = weights[ipt] * vrho[ipt];

	      const double z = zmat[mu_i]; // Z = N * B

	      const double dbx = dbasis_x_eval[mu_i]; // B_x
	      const double dby = dbasis_y_eval[mu_i]; // B_y
	      const double dbz = dbasis_z_eval[mu_i]; // B_z

	      g_acc_x += vrho_ipt * z * dbx;
	      g_acc_y += vrho_ipt * z * dby;
	      g_acc_z += vrho_ipt * z * dbz;

	      if( func.is_gga() or func.is_mgga() ) {
      	  // GGA Contributions
          const double vgamma_ipt = weights[ipt] * vgamma[ipt];

          const double dden_x = dden_x_eval[ipt];
          const double dden_y = dden_y_eval[ipt];
          const double dden_z = dden_z_eval[ipt];

	        const double zx = zmat_x[mu_i]; // Z_x = N * B_x
	        const double zy = zmat_y[mu_i]; // Z_y = N * B_y
	        const double zz = zmat_z[mu_i]; // Z_z = N * B_z

	        const double d2bxx = d2basis_xx_eval[mu_i]; // B^2_xx
	        const double d2bxy = d2basis_xy_eval[mu_i]; // B^2_xy
	        const double d2bxz = d2basis_xz_eval[mu_i]; // B^2_xz
	        const double d2byy = d2basis_yy_eval[mu_i]; // B^2_yy
	        const double d2byz = d2basis_yz_eval[mu_i]; // B^2_yz
	        const double d2bzz = d2basis_zz_eval[mu_i]; // B^2_zz

	        // sum_j B^2_{ij} * d_j n
	        double d2_term_x = d2bxx * dden_x + d2bxy * dden_y + d2bxz * dden_z;
	        double d2_term_y = d2bxy * dden_x + d2byy * dden_y + d2byz * dden_z;
	        double d2_term_z = d2bxz * dden_x + d2byz * dden_y + d2bzz * dden_z;

	        // sum_j (d_j n) * Z^j
	        double d11_zmat_term = dden_x * zx + dden_y * zy + dden_z * zz;

	        g_acc_x += 2 * vgamma_ipt * ( z * d2_term_x + dbx * d11_zmat_term );
	        g_acc_y += 2 * vgamma_ipt * ( z * d2_term_y + dby * d11_zmat_term );
	        g_acc_z += 2 * vgamma_ipt * ( z * d2_term_z + dbz * d11_zmat_term );
	      }
#if 0
	      if( func.is_mgga() ) {

                const double vtau_ipt = 0.5 * weights[ipt] * vtau[ipt];
	        const double zx = zmat_x[mu_i]; // Z_x = N * B_x
	        const double zy = zmat_y[mu_i]; // Z_y = N * B_y
	        const double zz = zmat_z[mu_i]; // Z_z = N * B_z
	        const double d2bxx = d2basis_xx_eval[mu_i]; // B^2_xx
	        const double d2bxy = d2basis_xy_eval[mu_i]; // B^2_xy
	        const double d2bxz = d2basis_xz_eval[mu_i]; // B^2_xz
	        const double d2byy = d2basis_yy_eval[mu_i]; // B^2_yy
	        const double d2byz = d2basis_yz_eval[mu_i]; // B^2_yz
	        const double d2bzz = d2basis_zz_eval[mu_i]; // B^2_zz
		double d2_term_x = d2bxx * zx + d2bxy * zy + d2bxz * zz;
		double d2_term_y = d2bxy * zx + d2byy * zy + d2byz * zz;
		double d2_term_z = d2bxz * zx + d2byz * zy + d2bzz * zz;

		g_acc_x += vtau_ipt * d2_term_x;
		g_acc_y += vtau_ipt * d2_term_y;
		g_acc_z += vtau_ipt * d2_term_z;

		if ( needs_laplacian ) {
		  const double vlapl_ipt = weights[ipt] * vlapl[ipt];
		  const double lbf = lbasis_eval[mu_i];
                  const double dlbx = dlbasis_x_eval[mu_i];
                  const double dlby = dlbasis_y_eval[mu_i];
                  const double dlbz = dlbasis_z_eval[mu_i];
		  d2_term_x = z * dlbx + zx * lbf + 2.0*d2_term_x;
		  d2_term_y = z * dlby + zy * lbf + 2.0*d2_term_y;
		  d2_term_z = z * dlbz + zz * lbf + 2.0*d2_term_z;

		  g_acc_x += vlapl_ipt * d2_term_x;
		  g_acc_y += vlapl_ipt * d2_term_y;
		  g_acc_z += vlapl_ipt * d2_term_z;

		}

	      }
#endif

      } // loop over bfns + grid points

      #pragma omp atomic
      EXC_GRAD[3*iAt + 0] += -2 * g_acc_x;
      #pragma omp atomic
      EXC_GRAD[3*iAt + 1] += -2 * g_acc_y;
      #pragma omp atomic
      EXC_GRAD[3*iAt + 2] += -2 * g_acc_z;

      bf_off += sh_sz; // Increment basis offset

    } // End loop over shells 
        
  } // End loop over tasks

  } // OpenMP Region

  
}

} // namespace GauXC::detail
