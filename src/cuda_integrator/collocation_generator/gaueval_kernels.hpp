#include <gauxc/shell.hpp>

namespace GauXC {

__global__
void gaueval_device_kernel(
  size_t         nshells,
  size_t         nbf,
  size_t         npts,
  const Shell*   shells_device,
  const size_t*  offs_device,
  const double*  pts_device,
  double*        eval_device
); 

__global__
void gaueval_device_kernel_deriv1(
  size_t         nshells,
  size_t         nbf,
  size_t         npts,
  const Shell*   shells_device,
  const size_t*  offs_device,
  const double*  pts_device,
  double*        eval_device,
  double*        deval_device_x,
  double*        deval_device_y,
  double*        deval_device_z
);

}
