#include "device/common/xc_functional_eval_wrapper.hpp"
#include "device_specific/cuda_util.hpp"

namespace GauXC {

void eval_kern_exc_vxc_lda( const functional_type& func, size_t npts,
  const double* rho, double* eps, double* vrho, type_erased_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  func.eval_exc_vxc_device( npts, rho, eps, vrho, stream );

}

void eval_kern_exc_vxc_gga( const functional_type& func, size_t npts,
  const double* rho, const double* gamma, double* eps, double* vrho,
  double* vgamma, type_erased_queue queue ) {

  cudaStream_t stream = queue.queue_as<util::cuda_stream>();
  func.eval_exc_vxc_device( npts, rho, gamma, eps, vrho, vgamma, stream );

}

}