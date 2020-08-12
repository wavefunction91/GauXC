#include <gauxc/shell.hpp>

namespace GauXC      {
namespace integrator {
namespace cuda       {

template <typename T>
void gaueval_device(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  cudaStream_t    stream
);

template <typename T>
void gaueval_device_deriv1(
  size_t          nshells,
  size_t          nbf,
  size_t          npts,
  const Shell<T>* shells_device,
  const size_t*   offs_device,
  const T*        pts_device,
  T*              eval_device,
  T*              deval_device_x,
  T*              deval_device_y,
  T*              deval_device_z,
  cudaStream_t    stream
);

} // namespace cuda
} // namespace integrator
} // namespace GauXC
