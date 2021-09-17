#pragma once
#include <gauxc/types.hpp>
#include "device/type_erased_queue.hpp"

namespace GauXC {

void eval_kern_exc_vxc_lda( const functional_type& func, size_t npts,
  const double* rho, double* eps, double* vrho, type_erased_queue queue );
void eval_kern_exc_vxc_gga( const functional_type& func, size_t npts,
  const double* rho, const double* gamma, double* eps, double* vrho,
  double* vgamma, type_erased_queue queue );

}