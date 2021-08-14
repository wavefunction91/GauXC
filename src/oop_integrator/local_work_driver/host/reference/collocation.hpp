#pragma once

#include <gauxc/basisset.hpp>

namespace GauXC {

void gau2grid_collocation( size_t                  npts, 
                           size_t                  nshells,
                           size_t                  nbe,
                           const double*           points, 
                           const BasisSet<double>& basis,
                           const int32_t*          shell_mask,
                           double*                 basis_eval );

void gau2grid_collocation_gradient( size_t                  npts, 
                                    size_t                  nshells,
                                    size_t                  nbe,
                                    const double*           points, 
                                    const BasisSet<double>& basis,
                                    const int32_t*          shell_mask,
                                    double*                 basis_eval, 
                                    double*                 dbasis_x_eval, 
                                    double*                 dbasis_y_eval,
                                    double*                 dbasis_z_eval );


}
