#pragma once

#include <gauxc/basisset.hpp>

namespace GauXC::integrator::host {

void eval_collocation( size_t                  npts, 
                       size_t                  nshells,
                       size_t                  nbe,
                       const double*           points, 
                       const BasisSet<double>& basis,
                       const int32_t*          shell_mask,
                       double*                 basis_eval );

void eval_collocation_deriv1( size_t                  npts, 
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
