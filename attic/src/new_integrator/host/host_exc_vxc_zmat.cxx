#include "host/host_exc_vxc_zmat.hpp"
#include "host/blas.hpp"

namespace GauXC  {
namespace integrator::host {

template <typename F>
void zmat_lda_host( int32_t   npts,
                    int32_t   nbf,
                    const F*  vrho,
                    const F*  basis,
                    F*        z_matrix ) {

  GauXC::blas::lacpy( 'A', nbf, npts, basis, nbf, 
                      z_matrix, nbf );

  for( int32_t i = 0; i < npts; ++i ) {

    auto* z_col = z_matrix + i*nbf;

    const F fact = 0.5 * vrho[i];
    GauXC::blas::scal( nbf, fact, z_col, 1 );

  }

} 

template
void zmat_lda_host( int32_t    npts,
                    int32_t    nbf,
                    const float*  vrho,
                    const float*  basis,
                    float*        z_matrix ); 
template
void zmat_lda_host( int32_t    npts,
                    int32_t    nbf,
                    const double*  vrho,
                    const double*  basis,
                    double*        z_matrix ); 



template <typename F>
void zmat_gga_host( int32_t   npts,
                    int32_t   nbf,
                    const F*  vrho,
                    const F*  vgamma,
                    const F*  basis,
                    const F*  dbasis_x,
                    const F*  dbasis_y,
                    const F*  dbasis_z,
                    const F*  dden_x,
                    const F*  dden_y,
                    const F*  dden_z,
                    F*        z_matrix ) {

  GauXC::blas::lacpy( 'A', nbf, npts, basis, nbf, 
                      z_matrix, nbf );

  for( int32_t i = 0; i < npts; ++i ) {

    const int32_t ioff = i * nbf;

    auto* z_col    = z_matrix + ioff;
    auto* bf_x_col = dbasis_x + ioff; 
    auto* bf_y_col = dbasis_y + ioff; 
    auto* bf_z_col = dbasis_z + ioff; 

    const F lda_fact = 0.5 * vrho[i];
    GauXC::blas::scal( nbf, lda_fact, z_col, 1 );

    const F gga_fact = 2. * vgamma[i]; 
    const auto x_fact = gga_fact * dden_x[i];
    const auto y_fact = gga_fact * dden_y[i];
    const auto z_fact = gga_fact * dden_z[i];

    GauXC::blas::axpy( nbf, x_fact, bf_x_col, 1, z_col, 1 );
    GauXC::blas::axpy( nbf, y_fact, bf_y_col, 1, z_col, 1 );
    GauXC::blas::axpy( nbf, z_fact, bf_z_col, 1, z_col, 1 );

  }

} 

template 
void zmat_gga_host( int32_t    npts,
                    int32_t    nbf,
                    const float*  vrho,
                    const float*  vgamma,
                    const float*  basis,
                    const float*  dbasis_x,
                    const float*  dbasis_y,
                    const float*  dbasis_z,
                    const float*  dden_x,
                    const float*  dden_y,
                    const float*  dden_z,
                    float*        z_matrix );

template 
void zmat_gga_host( int32_t    npts,
                    int32_t    nbf,
                    const double*  vrho,
                    const double*  vgamma,
                    const double*  basis,
                    const double*  dbasis_x,
                    const double*  dbasis_y,
                    const double*  dbasis_z,
                    const double*  dden_x,
                    const double*  dden_y,
                    const double*  dden_z,
                    double*        z_matrix );

}
}

