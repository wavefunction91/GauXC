
#include "host/reference_local_host_work_driver.hpp"
#include "host/reference/weights.hpp"
#include "host/reference/collocation.hpp"

#include "host/util.hpp"
#include "host/blas.hpp"
#include <stdexcept>

#include <gauxc/basisset_map.hpp>
#include "rys_integral.h"
#include "obara_saika_integrals.hpp"
#include "chebyshev_boys_computation.hpp"
#include <gauxc/util/real_solid_harmonics.hpp>
#include "integrator_util/integral_bounds.hpp"

namespace GauXC {

  ReferenceLocalHostWorkDriver::ReferenceLocalHostWorkDriver() {
    this -> boys_table = XCPU::boys_init();
  };
  
  ReferenceLocalHostWorkDriver::~ReferenceLocalHostWorkDriver() noexcept {
    XCPU::boys_finalize(this -> boys_table);
  };

  // Partition weights
  void ReferenceLocalHostWorkDriver::partition_weights( XCWeightAlg weight_alg, 
							const Molecule& mol, const MolMeta& meta, task_iterator task_begin, 
							task_iterator task_end ) {
    switch( weight_alg ) {
    case XCWeightAlg::Becke:
      reference_becke_weights_host( mol, meta, task_begin, task_end );
      break;
    case XCWeightAlg::SSF:
      reference_ssf_weights_host( mol, meta, task_begin, task_end );
      break;
    case XCWeightAlg::LKO:
      reference_lko_weights_host( mol, meta, task_begin, task_end );
      break;
    default:
      GAUXC_GENERIC_EXCEPTION("Weight Alg Not Supported");
    }
  }


  // Collocation
  void ReferenceLocalHostWorkDriver::eval_collocation( size_t npts, size_t nshells, 
						       size_t nbe, const double* pts, const BasisSet<double>& basis, 
						       const int32_t* shell_list, double* basis_eval ) {
    gau2grid_collocation( npts, nshells, nbe, pts, basis, shell_list, basis_eval );
  }


  // Collocation Gradient
  void ReferenceLocalHostWorkDriver::eval_collocation_gradient( size_t npts, 
								size_t nshells, size_t nbe, const double* pts, const BasisSet<double>& basis, 
								const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
								double* dbasis_y_eval, double* dbasis_z_eval) {
    gau2grid_collocation_gradient(npts, nshells, nbe, pts, basis, shell_list,
				  basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval );
  }

  void ReferenceLocalHostWorkDriver::eval_collocation_hessian( size_t npts, 
							       size_t nshells, size_t nbe, const double* pts, const BasisSet<double>& basis, 
							       const int32_t* shell_list, double* basis_eval, double* dbasis_x_eval, 
							       double* dbasis_y_eval, double* dbasis_z_eval, double* d2basis_xx_eval, 
							       double* d2basis_xy_eval, double* d2basis_xz_eval, double* d2basis_yy_eval, 
							       double* d2basis_yz_eval, double* d2basis_zz_eval ) {
    gau2grid_collocation_hessian(npts, nshells, nbe, pts, basis, shell_list,
				 basis_eval, dbasis_x_eval, dbasis_y_eval, dbasis_z_eval, d2basis_xx_eval,
				 d2basis_xy_eval, d2basis_xz_eval, d2basis_yy_eval, d2basis_yz_eval,
				 d2basis_zz_eval);
  }


  // X matrix (P * B)
  void ReferenceLocalHostWorkDriver::eval_xmat( size_t npts, size_t nbf, size_t nbe, 
						const submat_map_t& submat_map, const double* P, size_t ldp, 
						const double* basis_eval, size_t ldb, double* X, size_t ldx, double* scr ) {
    const auto* P_use = P;
    size_t ldp_use = ldp;
     
    if( submat_map.size() > 1 ) {
      detail::submat_set( nbf, nbf, nbe, nbe, P, ldp, scr, nbe, submat_map );
      P_use = scr;
      ldp_use = nbe;
    } else if( nbe != nbf ) {
      P_use = P + submat_map[0][0]*(ldp+1);
    }

    blas::gemm( 'N', 'N', nbe, npts, nbe, 2., P_use, ldp_use, basis_eval, ldb, 
		0., X, ldx );

  }

  // U/VVar LDA (density)
  void ReferenceLocalHostWorkDriver::eval_uvvar_lda( size_t npts, size_t nbe, 
						     const double* basis_eval, const double* X, size_t ldx, double* den_eval) {
    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioff = size_t(i) * ldx;
      const auto*   X_i = X + ioff;
      den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );
      
    }

  }

  //VVar GGA (density + grad, gamma)
  void ReferenceLocalHostWorkDriver::eval_uvvar_gga( size_t npts, size_t nbe, 
						     const double* basis_eval, const double* dbasis_x_eval, 
						     const double *dbasis_y_eval, const double* dbasis_z_eval, const double* X, 
						     size_t ldx, double* den_eval, double* dden_x_eval, double* dden_y_eval, 
						     double* dden_z_eval, double* gamma ) {

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const size_t ioff = size_t(i) * ldx;
      const auto*   X_i = X + ioff;

      den_eval[i] = blas::dot( nbe, basis_eval + ioff, 1, X_i, 1 );

      const auto dx = 
	2. * blas::dot( nbe, dbasis_x_eval + ioff, 1, X_i, 1 );
      const auto dy = 
	2. * blas::dot( nbe, dbasis_y_eval + ioff, 1, X_i, 1 );
      const auto dz = 
	2. * blas::dot( nbe, dbasis_z_eval + ioff, 1, X_i, 1 );

      dden_x_eval[i] = dx;
      dden_y_eval[i] = dy;
      dden_z_eval[i] = dz;

      gamma[i] = dx*dx + dy*dy + dz*dz;

    }
  }

  // Eval Z Matrix LDA VXC
  void ReferenceLocalHostWorkDriver::eval_zmat_lda_vxc( size_t npts, size_t nbf, 
							const double* vrho, const double* basis_eval, double* Z, size_t ldz ) {


    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Z, ldz );

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      auto* z_col = Z + i*ldz;

      const double fact = 0.5 * vrho[i];
      GauXC::blas::scal( nbf, fact, z_col, 1 );

    }

  }

  // Eval Z Matrix GGA VXC
  void ReferenceLocalHostWorkDriver::eval_zmat_gga_vxc( size_t npts, size_t nbf, 
							const double* vrho, const double* vgamma, const double* basis_eval, 
							const double* dbasis_x_eval, const double* dbasis_y_eval, 
							const double* dbasis_z_eval, const double* dden_x_eval, 
							const double* dden_y_eval, const double* dden_z_eval, double* Z, size_t ldz ) {

    if( ldz != nbf ) GAUXC_GENERIC_EXCEPTION(std::string("INVALID DIMS"));
    blas::lacpy( 'A', nbf, npts, basis_eval, nbf, Z, nbf );

    for( int32_t i = 0; i < (int32_t)npts; ++i ) {

      const int32_t ioff = i * nbf;

      auto* z_col    = Z + ioff;
      auto* bf_x_col = dbasis_x_eval + ioff; 
      auto* bf_y_col = dbasis_y_eval + ioff; 
      auto* bf_z_col = dbasis_z_eval + ioff; 

      const auto lda_fact = 0.5 * vrho[i];
      blas::scal( nbf, lda_fact, z_col, 1 );

      const auto gga_fact = 2. * vgamma[i]; 
      const auto x_fact = gga_fact * dden_x_eval[i];
      const auto y_fact = gga_fact * dden_y_eval[i];
      const auto z_fact = gga_fact * dden_z_eval[i];

      blas::axpy( nbf, x_fact, bf_x_col, 1, z_col, 1 );
      blas::axpy( nbf, y_fact, bf_y_col, 1, z_col, 1 );
      blas::axpy( nbf, z_fact, bf_z_col, 1, z_col, 1 );

    }

  }

  // Increment VXC by Z
  void ReferenceLocalHostWorkDriver::inc_vxc( size_t npts, size_t nbf, size_t nbe, 
					      const double* basis_eval, const submat_map_t& submat_map, const double* Z, 
					      size_t ldz, double* VXC, size_t ldvxc, double* scr ) {

    if( submat_map.size() > 1 ) {
      blas::syr2k('L', 'N', nbe, npts, 1., basis_eval, nbe, Z, ldz, 0., scr, nbe );
      detail::inc_by_submat( nbf, nbf, nbe, nbe, VXC, ldvxc, scr, nbe, submat_map );
    } else {
      blas::syr2k('L', 'N', nbe, npts, 1., basis_eval, nbe, Z, ldz, 1., 
		  VXC + submat_map[0][0]*(ldvxc+1), ldvxc );
    }

  }

  // Increment K by G
  void ReferenceLocalHostWorkDriver::inc_exx_k( size_t npts, size_t nbf, 
						size_t nbe_bra, size_t nbe_ket, const double* basis_eval, 
						const submat_map_t& submat_map_bra, const submat_map_t& submat_map_ket, 
						const double* G, size_t ldg, double* K, size_t ldk, double* scr ) {

    if( submat_map_bra.size() > 1 or submat_map_ket.size() > 1 ) {
      blas::gemm( 'N', 'T', nbe_bra, nbe_ket, npts, 1., basis_eval, nbe_bra,
		  G, ldg, 0., scr, nbe_bra );

      detail::inc_by_submat( nbf, nbf, nbe_bra, nbe_ket, K, ldk, scr, nbe_bra, 
			     submat_map_bra, submat_map_ket );
    } else {
      blas::gemm( 'N', 'T', nbe_bra, nbe_ket, npts, 1., basis_eval, nbe_bra,
		  G, ldg, 1., K + submat_map_ket[0][0]*ldk + submat_map_bra[0][0], ldk );
    }

  }

  struct RysBasis {
    std::vector< shells > _shells;
    RysBasis( const BasisSet<double>& basis ) {
      size_t nshells = basis.size();
      _shells.resize(nshells);
      for( int i = 0; i < nshells; ++i ) {
	_shells[i].origin.x = basis[i].O()[0];
	_shells[i].origin.y = basis[i].O()[1];
	_shells[i].origin.z = basis[i].O()[2];

	_shells[i].m = basis[i].nprim();
	_shells[i].L = basis[i].l();

	_shells[i].coeff = new coefficients[_shells[i].m];
	for( int j = 0; j < _shells[i].m; ++j ) {
	  _shells[i].coeff[j].alpha = basis[i].alpha()[j];
	  _shells[i].coeff[j].coeff = basis[i].coeff()[j];
	}
      }
    }

    shells& operator[](int i){ return _shells.at(i); }
    const shells& operator[](int i) const { return _shells.at(i); }

    ~RysBasis() noexcept {
      for( auto& sh : _shells ) delete sh.coeff;
    }
  };

  void ReferenceLocalHostWorkDriver::eval_exx_fmat( size_t npts, size_t nbf, 
						    size_t nbe_bra, size_t nbe_ket, const submat_map_t& submat_map_bra,
						    const submat_map_t& submat_map_ket, const double* P, size_t ldp,
						    const double* basis_eval, size_t ldb, double* F, size_t ldf,
						    double* scr ) {

    const auto* P_use = P;
    size_t ldp_use = ldp;

    if( submat_map_bra.size() > 1 or submat_map_ket.size() > 1 ) {
      detail::submat_set( nbf, nbf, nbe_bra, nbe_ket, P, ldp,
			  scr, nbe_bra, submat_map_bra, submat_map_ket );
      P_use = scr;
      ldp_use = nbe_bra;
    } else {
      P_use = P + submat_map_ket[0][0]*ldp + submat_map_bra[0][0];
    }

    blas::gemm( 'N', 'N', nbe_bra, npts, nbe_ket, 1., P_use, ldp_use, basis_eval,
		ldb, 0., F, ldf );

  }

  void ReferenceLocalHostWorkDriver:: eval_exx_gmat( size_t npts, size_t nshells,
						     size_t nbe, const double* points, const double* weights, 
						     const BasisSet<double>& basis, const BasisSetMap& basis_map, 
						     const int32_t* shell_list, const double* X, size_t ldx, double* G, size_t ldg ) {
    
    // Cast points to Rys format (binary compatable)
    point* _points = reinterpret_cast<point*>(const_cast<double*>(points));
    std::vector<double> _points_transposed(3 * npts);

    for(int i = 0; i < npts; ++i) {
      _points_transposed[i + 0 * npts] = _points[i].x;
      _points_transposed[i + 1 * npts] = _points[i].y;
      _points_transposed[i + 2 * npts] = _points[i].z;
    }
  
    // Set G to zero
    for( int j = 0; j < npts; ++j )
      for( int i = 0; i < nbe;  ++i ) {
	G[i + j*ldg] = 0.;
      }

    // Copy the basis set 
    RysBasis rys_basis(basis);

    size_t offset = 0;
    std::vector<shell_pair> shpairs;
    shpairs.resize(nshells * (nshells + 1) / 2);
    for( size_t i = 0; i < nshells; ++i ) {
      for( size_t j = 0; j <= i; ++j ) {
	if( rys_basis._shells[i].L >= rys_basis._shells[j].L )
	  XCPU::generate_shell_pair(rys_basis._shells[i], rys_basis._shells[j], shpairs[offset]);
	else
	  XCPU::generate_shell_pair(rys_basis._shells[j], rys_basis._shells[i], shpairs[offset]);

	offset++;
      }
    }

    // Spherical Harmonic Transformer
    util::SphericalHarmonicTransform sph_trans(5);

    const bool any_pure = std::any_of( shell_list, shell_list + nshells,
				       [&](const auto& i){ return basis.at(i).pure(); } );
    
    const size_t nbe_cart = basis.nbf_cart_subset( shell_list, shell_list + nshells );

    std::vector<double> X_cart, G_cart;
    if( any_pure ){
      X_cart.resize( nbe_cart * npts );
      G_cart.resize( nbe_cart * npts, 0. );

      // Transform X into cartesian
      int ioff = 0;
      int ioff_cart = 0;
      for( int i = 0; i < nshells; ++i ) {
	const auto ish = shell_list[i];
	const auto& shell      = basis.at(ish);
	const int shell_l       = shell.l();
	const int shell_sz      = shell.size();
	const int shell_cart_sz = shell.cart_size();

	if( shell.pure() and shell_l > 0 ) {
	  sph_trans.itform_bra_cm( shell_l, npts, X + ioff, ldx,
				   X_cart.data() + ioff_cart, nbe_cart );
	} else {
	  blas::lacpy( 'A', shell_sz, npts, X + ioff, ldx,
		       X_cart.data() + ioff_cart, nbe_cart );
	}
	ioff += shell_sz;
	ioff_cart += shell_cart_sz;
      }
    }

    const auto* X_use = any_pure ? X_cart.data() : X;
    auto*       G_use = any_pure ? G_cart.data() : G;
    const auto ldx_use = any_pure ? nbe_cart : ldx;
    const auto ldg_use = any_pure ? nbe_cart : ldg;

    std::vector<double> X_cart_rm( nbe_cart*npts,0. ), G_cart_rm( nbe_cart*npts,0. );
    for( auto i = 0; i < nbe_cart; ++i )
      for( auto j = 0; j < npts;     ++j ) {
	X_cart_rm[i*npts + j] = X_use[i + j*ldx_use];
      }

    offset = 0;
    size_t ioff_cart = 0;
    for( int i = 0; i < nshells; ++i ) {
      const auto ish        = shell_list[i];
      const auto& bra       = basis.at(ish);
      const int bra_cart_sz = bra.cart_size();

      size_t joff_cart = 0;
      for( int j = 0; j <= i; ++j ) {
	const auto jsh        = shell_list[j];
	const auto& ket       = basis.at(jsh);
	const int ket_cart_sz = ket.cart_size();

	XCPU::compute_integral_shell_pair( npts, ish == jsh, bra.L, ket.L, shpairs[offset], _points_transposed.data(),
					   X_cart_rm.data()+ioff_cart*npts, X_cart_rm.data()+joff_cart*npts, npts,
					   G_cart_rm.data()+ioff_cart*npts, G_cart_rm.data()+joff_cart*npts, npts,
					   const_cast<double*>(weights), this -> boys_table );

	offset++;
	joff_cart += ket_cart_sz;
      }
	
      ioff_cart += bra_cart_sz;
    }
   
    for( auto i = 0; i < nbe_cart; ++i )
      for( auto j = 0; j < npts;     ++j ) {
	G_use[i + j*ldx_use] = G_cart_rm[i*npts + j];
      }
  
    // Transform G back to spherical
    if( any_pure ) {
      int ioff = 0;
      int ioff_cart = 0;
      for( int i = 0; i < nshells; ++i ) {
	const auto ish = shell_list[i];
	const auto& shell      = basis.at(ish);
	const int shell_l       = shell.l();
	const int shell_sz      = shell.size();
	const int shell_cart_sz = shell.cart_size();

	if( shell.pure() and shell_l > 0 ) {
	  sph_trans.tform_bra_cm( shell_l, npts, G_cart.data() + ioff_cart, nbe_cart,
				  G + ioff, ldg );
	} else {
	  blas::lacpy( 'A', shell_sz, npts, G_cart.data() + ioff_cart, nbe_cart,
		       G + ioff, ldg );
	}
	ioff += shell_sz;
	ioff_cart += shell_cart_sz;
      }
    }
  }

}
