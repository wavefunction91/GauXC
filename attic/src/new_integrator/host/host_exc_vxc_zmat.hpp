#pragma once
#include <cstdint>

namespace GauXC  {
namespace integrator::host {

template <typename F>
void zmat_lda_host( int32_t   npts,
                    int32_t   nbf,
                    const F*  vrho,
                    const F*  basis,
                    F*        z_matrix ); 

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
                    F*        z_matrix ); 

}
}
