#include "host_collocation.hpp"
#include "gau2grid.h"

namespace GauXC::integrator::host {

void eval_collocation( size_t                  npts, 
                       size_t                  nshells,
                       size_t                  nbe,
                       const double*           points, 
                       const BasisSet<double>& basis,
                       const int32_t*          shell_mask,
                       double*                 basis_eval ) {

  std::allocator<double> a;
  auto* rv = a.allocate( npts * nbe );

  size_t ncomp = 0;
  for( int i = 0; i < nshells; ++i ) {

    const auto& sh = basis.at(shell_mask[i]);
    int order = sh.pure() ? GG_SPHERICAL_CCA : GG_CARTESIAN_CCA; 
    gg_collocation( sh.l(), npts, points, 3, sh.nprim(), sh.coeff_data(),
      sh.alpha_data(), sh.O_data(), order, rv + ncomp*npts );

    ncomp += sh.size();

  }

  gg_fast_transpose( ncomp, npts, rv, basis_eval );
  a.deallocate( rv, npts*nbe );

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

  std::allocator<double> a;
  auto* rv = a.allocate( 4 * npts * nbe );
  auto* rv_x = rv   + npts * nbe;
  auto* rv_y = rv_x + npts * nbe;
  auto* rv_z = rv_y + npts * nbe;

  size_t ncomp = 0;
  for( int i = 0; i < nshells; ++i ) {

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

}


}
