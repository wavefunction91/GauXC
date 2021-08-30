#include "host/host_collocation.hpp"


#ifdef GAUXC_ENABLE_GAU2GRID
  #include "gau2grid/gau2grid.h"
#else
  #include "collocation/collocation_angular_cartesian.hpp"
  #include "collocation/collocation_angular_spherical_unnorm.hpp"
  #include "collocation/collocation_radial.hpp"
#endif

namespace GauXC::integrator::host {

void eval_collocation( size_t                  npts, 
                       size_t                  nshells,
                       size_t                  nbe,
                       const double*           points, 
                       const BasisSet<double>& basis,
                       const int32_t*          shell_mask,
                       double*                 basis_eval ) {

#ifdef GAUXC_ENABLE_GAU2GRID

  std::allocator<double> a;
  auto* rv = a.allocate( npts * nbe );

  size_t ncomp = 0;
  for( size_t i = 0; i < nshells; ++i ) {

    const auto& sh = basis.at(shell_mask[i]);
    int order = sh.pure() ? GG_SPHERICAL_CCA : GG_CARTESIAN_CCA; 
    gg_collocation( sh.l(), npts, points, 3, sh.nprim(), sh.coeff_data(),
      sh.alpha_data(), sh.O_data(), order, rv + ncomp*npts );

    ncomp += sh.size();

  }

  gg_fast_transpose( ncomp, npts, rv, basis_eval );
  a.deallocate( rv, npts*nbe );

#else
  
  for( size_t ipt = 0; ipt < npts;  ++ipt )
  for( size_t i = 0;   i < nshells; ++i   ) {
    
    const auto ish = shell_mask[i];
    const auto& sh = basis.at(ish);
    auto* eval = basis_eval + ipt*nbe + basis.shell_to_first_ao( ish );

    double x,y,z, bf;
    integrator::cuda::collocation_device_radial_eval( sh, points + 3*ipt, 
                                                      &x, &y, &z, &bf );

    if( sh.pure() )
      integrator::cuda::collocation_spherical_unnorm_angular( sh.l(), bf, x, y, z,
                                                              eval );
    else
      integrator::cuda::collocation_cartesian_angular( sh.l(), bf, x, y, z, eval );
                                                              
                                                              
  }

#endif

}

void eval_collocation_deriv1( size_t                  npts, 
                              size_t                  nshells,
                              size_t                  nbe,
                              const double*           points, 
                              const BasisSet<double>& basis,
                              const int32_t*          shell_mask,
                              double*                 basis_eval, 
                              double*                 dbasis_x_eval, 
                              double*                 dbasis_y_eval,
                              double*                 dbasis_z_eval ) {

#ifdef GAUXC_ENABLE_GAU2GRID

  std::allocator<double> a;
  auto* rv = a.allocate( 4 * npts * nbe );
  auto* rv_x = rv   + npts * nbe;
  auto* rv_y = rv_x + npts * nbe;
  auto* rv_z = rv_y + npts * nbe;

  size_t ncomp = 0;
  for( size_t i = 0; i < nshells; ++i ) {

    const auto& sh = basis.at(shell_mask[i]);
    int order = sh.pure() ? GG_SPHERICAL_CCA : GG_CARTESIAN_CCA; 
    gg_collocation_deriv1( sh.l(), npts, points, 3, sh.nprim(), sh.coeff_data(),
      sh.alpha_data(), sh.O_data(), order, rv + ncomp*npts, 
      rv_x + ncomp*npts, rv_y + ncomp*npts, rv_z + ncomp*npts );

    ncomp += sh.size();

  }

  gg_fast_transpose( ncomp, npts, rv,   basis_eval );
  gg_fast_transpose( ncomp, npts, rv_x, dbasis_x_eval );
  gg_fast_transpose( ncomp, npts, rv_y, dbasis_y_eval );
  gg_fast_transpose( ncomp, npts, rv_z, dbasis_z_eval );

  a.deallocate( rv, 4*npts*nbe );

#else 

  for( size_t ipt = 0; ipt < npts;  ++ipt )
  for( size_t i = 0;   i < nshells; ++i   ) {
    
    const auto ish = shell_mask[i];
    const auto& sh = basis.at(ish);
    auto* eval = basis_eval + ipt*nbe + basis.shell_to_first_ao( ish );
    auto* deval_x = dbasis_x_eval + ipt*nbe + basis.shell_to_first_ao( ish );
    auto* deval_y = dbasis_y_eval + ipt*nbe + basis.shell_to_first_ao( ish );
    auto* deval_z = dbasis_z_eval + ipt*nbe + basis.shell_to_first_ao( ish );

    double x,y,z, bf, dbf_x, dbf_y, dbf_z;
    integrator::cuda::collocation_device_radial_eval_deriv1( sh, points + 3*ipt, 
                                                      &x, &y, &z, &bf, &dbf_x,
                                                      &dbf_y, &dbf_z);

    if( sh.pure() )
      integrator::cuda::collocation_spherical_unnorm_angular_deriv1( 
        sh.l(), bf, dbf_x, dbf_y, dbf_z, x, y, z, eval, deval_x, deval_y, deval_z );
    else
      integrator::cuda::collocation_cartesian_angular_deriv1( 
        sh.l(), bf, dbf_x, dbf_y, dbf_z, x, y, z, eval, deval_x, deval_y, deval_z );
                                                              
  }

#endif
}


}
